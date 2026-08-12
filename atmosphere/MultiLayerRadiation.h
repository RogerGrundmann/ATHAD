#pragma once

#include "MixtureAtm.h"
#include "cAtmosphereModel.h"
#include "Utils.h"

#include <vector>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

// ============================================================================
// MultiLayerRadiation — friend class of cAtmosphereModel
//
// Multi-layer long/short-wave radiation model for the surface-temperature
// distribution (computation of the local temperature from short- and long-wave
// radiation). For each (j,k) column it builds:
//   - the grey-body emission of every layer (Stefan-Boltzmann, sigma T^4),
//   - the layer emissivity epsilon from the water-vapour law (Bignami 1995,
//     valid -40..45 C) plus an optional CO2 term (Atwater & Ball) — the CO2
//     contribution is presently disabled (eps_co2 = 0; the CO2 coupling is
//     rough and to be updated),
// then solves the layer radiative-balance as a tridiagonal system with the
// Thomas algorithm and inverts the resulting radiation back to temperature.
// The surface layer receives the albedo-reduced incoming short-wave flux.
//
// Refactored from the former free member function
// cAtmosphereModel::RadiationMultiLayer() into the project's friend-class idiom
// (mirrors ThermoAtm / TurbulenceAtm / PressureSolverAtm). The (j,k) columns are
// independent, so the outer j loop is OpenMP-parallel with ALL Thomas-solver
// scratch (step/alfa/beta/AA/CA/CC and the per-column radiation_original)
// thread-local — the original reused a single shared set of scratch vectors,
// which would race under OpenMP.
//
// Reads : t, c, co2, p_stat, h, layer heights, albedo/short-wave/emissivity
//         constants (albedo_pole/equator, rad_pole/equator_short, sigma, ep,
//         co2_0, t_0, p_0).
// Writes: t (updated), radiation, epsilon, epsilon_2D, albedo, short_wave_radiation.
//
// Usage:  MultiLayerRadiation(*this).run();
// ============================================================================


class MultiLayerRadiation {
public:
    explicit MultiLayerRadiation(cAtmosphereModel& model)
        : m(model)
    {}

    void run()
    {
        using namespace std;
        cout << endl << "      RadiationMultiLayer" << endl;
        auto begin = std::chrono::high_resolution_clock::now();

        const int j_max  = m.jm - 1;
        const int j_half = j_max / 2;

        // ---- latitude profiles (computed once; read-only in the parallel loop) ----

        // Surface albedo with a temperature-dependent ICE/SNOW FEEDBACK (replaces the former
        // fixed pole->equator parabola, which froze the feedback out). Base value by surface
        // type (ocean vs land); where the surface is cold, blend toward a bright snow/sea-ice
        // albedo over a smooth ramp around freezing. Recomputed every MLR call on the CURRENT
        // surface T, so it is a LIVE feedback: warming -> less ice -> lower albedo -> more
        // absorbed shortwave -> extra (polar-amplified) warming that the dynamics cannot mix
        // away (unlike a forcing/nudge). Constants are tunable. NOTE the cloud SW bump below
        // overwrites this per column where cloud is present, so its net reach is cloud-limited.
        // ATHAD: the surface is a magma ocean, and there is no ice-albedo feedback.
        //
        // What was here: an ocean/land base albedo with a ramp toward snow/sea-ice as the
        // surface cooled through 275 -> 265 K. None of it applies. There is no land, no
        // snow, and at 1500 K nothing within 1200 K of the ice thresholds — the ramp was
        // dead code that always returned the ice-free ocean value.
        //
        // A quenching silicate melt is dark: measured basaltic-melt albedos are ~0.05-0.10.
        // This is the CLEAR-SKY value only; the cloud bump below raises it wherever the
        // model actually produces condensate. That separation is the point of this fix —
        // the reflective cloud deck a runaway greenhouse is supposed to have must be EARNED
        // by condensate the model generated, not asserted as a constant albedo. Asserting
        // 0.4 while the column condenses nothing was the inconsistency being removed.
        //
        // Now a PARAMETER (albedo_surface), not a literal. It was written here as a
        // constexpr while the config carried an inert albedo_pole/albedo_equator pair that
        // nothing read — so the file said one thing and the configuration said another.
        const double alb_surface_molten = m.albedo_surface;   // dark silicate melt, clear sky
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.albedo.y[j][k] = alb_surface_molten;

        // Incoming short-wave radiation: pole -> equator parabola, hemispherically symmetric.
        m.short_wave_radiation = std::vector<double>(m.jm, m.rad_pole_short);
        const double rad_short_eff = m.rad_pole_short - m.rad_equator_short;
        for (int j = j_half; j >= 0; j--)
            m.short_wave_radiation[j] =
                rad_short_eff * parabola((double)j / (double)j_half) + m.rad_pole_short;
        for (int j = j_max; j > j_half; j--)
            m.short_wave_radiation[j] = m.short_wave_radiation[j_max - j];

        // ---- per-column radiative balance (columns independent -> OpenMP over j) ----
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < m.jm; j++) {

            // Thread-local Thomas-solver scratch; reused across k within this j.
            // Every column fully overwrites the entries it later reads (i_mount = 0,
            // i_trop = im-1), so reuse is race-free.
            std::vector<double> alfa(m.im, 0.0), beta(m.im, 0.0);
            std::vector<double> AA(m.im, 0.0), CA(m.im, 0.0);
            std::vector<double> radiation_original(m.im, 0.0);
            std::vector<std::vector<double> > CC(m.im, std::vector<double>(m.im, 0.0));

            const int i_trop  = m.im - 1;   // top layer (tropopause proxy)
            const int i_mount = 0;          // surface / bottom layer

            for (int k = 0; k < m.km; k++) {

                // Grey-body emission of each layer and its "original" reference.
                for (int i = i_mount; i <= i_trop; i++) {
                    m.radiation.x[i][j][k] = m.sigma * pow(m.t.x[i][j][k] * m.t_0, 4.0);
                    radiation_original[i]  = m.radiation.x[i][j][k];
                }

                // ATHAD: layer optical depth from COLUMN MASS with pressure broadening.
                //
                // What was here: the Bignami (1995) clear-sky column emissivity
                // eps = 0.684 + 0.0056*e_surf, split into a "dry baseline" and a
                // water-vapour part and distributed over the layers, plus an Atwater & Ball
                // CO2 band with a scale factor tuned to 0.17 to reproduce Earth's ~30 W/m2
                // of CO2 greenhouse. Every one of those numbers is a regression on
                // present-day terrestrial columns — Bignami is a fit to MEDITERRANEAN SEA
                // SURFACE measurements. At 250 bar with a 67 %-by-mass water column they are
                // extrapolated some four orders of magnitude beyond their calibration, and
                // eps = 0.684 + 0.0056*e_surf saturates to 0.999 the instant e_surf exceeds
                // ~56 hPa (here it is ~2e5 hPa), so the whole scheme degenerates to "every
                // layer is a blackbody" and carries no information about the composition.
                //
                // What replaces it: a grey optical depth built from the absorber mass each
                // layer actually contains,
                //
                //     tau_i = SUM_s kappa_s * u_s,i * (p_i / p_ref)
                //     u_s,i = q_s * dp_i / g            [kg/m2]  column mass of species s
                //
                // The (p_i/p_ref) factor is pressure broadening: collisional line widths grow
                // in proportion to pressure, so the absorption per unit mass does too. It is
                // the term that matters most here and the one no Earth-calibrated emissivity
                // fit contains — at 250 bar it is a factor of 250 over the 1 bar reference,
                // and it is what makes the deep atmosphere opaque and the thin top
                // transparent, rather than a fit saturating everywhere.
                //
                // kappa values are grey-band mass absorption coefficients [m2/kg]. The water
                // value is the one conventionally used for grey runaway-greenhouse models;
                // CO2 is weaker per unit mass; the N2/CO/CH4 background is nearly
                // transparent in the thermal infrared. These are ASSUMPTIONS with roughly a
                // factor-of-two uncertainty and are the single biggest lever on the answer —
                // see README, and the OLR check below, which is what tests them.
                const double kappa_H2O = m.kappa_H2O;                 // [m2/kg]
                const double kappa_CO2 = m.kappa_CO2;                 // [m2/kg]
                const double kappa_bg  = m.kappa_bg;                  // [m2/kg]
                constexpr double p_ref = 1.0e5;                       // [Pa] 1 bar broadening reference
                const double inv_g     = 1.0 / m.g;

                double lwp_col = 0.0, iwp_col = 0.0;                  // condensate paths [g/m2] (SW albedo bump)
                for (int i = i_mount; i <= i_trop; i++) {
                    // Layer mass, from the pressure drop across it. dp in hPa -> Pa.
                    double dp_Pa = (i < i_trop)
                                 ? (m.p_stat.x[i][j][k] - m.p_stat.x[i+1][j][k]) * 100.0
                                 :  m.p_stat.x[i][j][k] * 100.0;      // top layer carries all mass above
                    if (dp_Pa < 0.0) dp_Pa = 0.0;

                    const double q_v = std::max(0.0, m.c.x[i][j][k]);
                    const double q_c = std::max(0.0, m.co2.x[i][j][k]);
                    const double q_b = std::max(0.0, 1.0 - q_v - q_c);

                    const double u_col = dp_Pa * inv_g;               // [kg/m2] total layer mass
                    const double broad = m.p_stat.x[i][j][k] * 100.0 / p_ref;   // pressure broadening

                    double tau_gas = (kappa_H2O * q_v + kappa_CO2 * q_c + kappa_bg * q_b)
                                   * u_col * broad;

                    // Cloud liquid + ice longwave greenhouse, unchanged in form (Stephens
                    // 1978 mass absorption), but the density now comes from the LOCAL mixture
                    // gas constant instead of the hard-coded 287.0 J/(kg K) of dry Earth air —
                    // which is 26 % off here and was applied to every layer.
                    constexpr double k_liq = 0.12, k_ice = 0.055;     // LW mass absorption [m2/g]
                    const double dz  = (i < i_trop) ? (m.get_layer_height(i+1) - m.get_layer_height(i))
                                                    : (m.get_layer_height(i) - m.get_layer_height(i-1));
                    const double T_i   = m.t.x[i][j][k] * m.t_0;
                    const double R_i   = AtmMixture::R_of(q_v, q_c, m.m_comp.R_bg);
                    const double rho_i = (T_i > 0.0) ? (m.p_stat.x[i][j][k] * 100.0) / (R_i * T_i) : 0.0;
                    const double cw_l  = std::max(0.0, m.cloud.x[i][j][k]);
                    const double cw_i  = std::max(0.0, m.ice.x[i][j][k]);
                    const double LWP_i = cw_l * rho_i * dz * 1000.0;  // [g/m2]
                    const double IWP_i = cw_i * rho_i * dz * 1000.0;  // [g/m2]
                    lwp_col += LWP_i;  iwp_col += IWP_i;

                    const double tau = tau_gas + k_liq * LWP_i + k_ice * IWP_i;

                    // exp(-tau) underflows to 0 for the huge optical depths the deep column
                    // carries, which is the correct answer (eps = 1, a perfect blackbody
                    // layer) — but guard it so no NaN can come out of the exponential.
                    m.epsilon.x[i][j][k] = (tau > 700.0) ? 1.0 : (1.0 - exp(-tau));
                }
                m.epsilon_2D.y[j][k] = m.epsilon.x[i_mount][j][k];

                // Cloud/ice SHORTWAVE albedo bump (stage 2): reflective clouds raise the column
                // albedo toward a cloud value, cutting absorbed SW (cooling) to compete with the LW
                // greenhouse above. alpha_eff = alpha_surf + (alpha_cloud - alpha_surf)*
                // (1 - exp(-k_sw*CWP_sw)), CWP_sw = LWP + f_ice*IWP (ice is optically thinner / less
                // reflective per unit path, so weighted down -> thin cirrus stays net-warming while
                // thick low liquid cloud net-cools). Overwrites the clear-sky latitude albedo for
                // this column; it is read below by the SW terms (tridiagonal source dd + surface
                // energy balance) and re-derived fresh each MLR call. See project_multilayer_radiation.
                {
                    // Two-stream-like cloud SW reflectivity on a CAPPED condensate path. The old
                    // exp(-k_sw*CWP) with k_sw=0.030 saturated at CWP~100 g/m2, so — and especially
                    // with the model's excessive cloud water (LWP ~1500 g/m2 vs observed ~100) —
                    // EVERY cloudy column was pinned at the asymptotic 0.60, wiping out the latitude
                    // gradient and MASKING the surface ice-albedo feedback (poles read the same 0.60
                    // as tropical ocean). Instead cap the path the radiation sees at a physical
                    // thick-cloud value and let reflectivity rise GENTLY as refl = tau/(tau+2)
                    // (0.5 at tau=2), composited over the surface (ice-feedback) albedo. Cloudy
                    // tropics now land ~0.3, thin/clear cells relax toward the surface value, and
                    // the polar ice albedo shows through — a physical gradient. (The same cloud-water
                    // excess still inflates the LW tau_cloud above; capping that is the next step.)
                    // alpha_cloud is a PARAMETER (albedo_cloud), and on ATHAD it is very
                    // nearly the entire planetary albedo: the deck's condensate path runs to
                    // 1e5 g/m2 against a cwp_tau of 100, so refl = tau/(tau+2) saturates and
                    // every cloudy column returns alpha_cloud to four decimals. Whatever
                    // this number is set to IS the model's albedo. It was a bare literal.
                    const double alpha_cloud     = m.albedo_cloud;   // thick cloud-top SW albedo
                    constexpr double f_ice_sw    = 0.50;   // ice SW reflectivity weight vs liquid
                    constexpr double cwp_tau     = 100.0;  // g/m2 per unit effective optical thickness
                    const double cwp_sw = lwp_col + f_ice_sw * iwp_col;          // [g/m2]
                    const double tau    = cwp_sw / cwp_tau;
                    const double refl   = tau / (tau + 2.0);                     // gentle saturation (0.5 at tau=2)
                    const double a0     = m.albedo.y[j][k];                      // surface (ice-feedback) albedo
                    if (alpha_cloud > a0)
                        m.albedo.y[j][k] = a0 + (alpha_cloud - a0) * refl;
                }

                // Transmitted (AA) / absorbed (CC diagonal) radiation, and the sum CA
                // of all radiations transmitted through each layer.
                AA[i_mount]          = m.radiation.x[i_mount][j][k];              // surface radiation
                CC[i_mount][i_mount] = m.epsilon.x[i_mount][j][k] * m.radiation.x[i_mount][j][k];
                for (int i = i_mount + 1; i <= i_trop; i++) {
                    AA[i]    = AA[i - 1] * (1.0 - m.epsilon.x[i][j][k]);          // transmitted from each layer
                    CC[i][i] = m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];     // absorbed in each layer
                }
                for (int i = i_mount + 2; i <= i_trop; i++) {
                    CA[i] = 0.0;
                    for (int l = 1; l <= i - 1; l++) {
                        CC[l][i] = CC[l][i - 1] * (1.0 - m.epsilon.x[i][j][k]);   // transmitted past layer i
                        CA[i] += CC[l][i];                                        // sum over all l
                    }
                }

                // Thomas algorithm — forward elimination (alfa/beta recurrence).
                // The sweep now INCLUDES the top row i = i_trop. The old code stopped at
                // i_trop-1 and closed the system with an ad-hoc top formula that divided by
                // (CA[i_trop]-CA[i_trop-1]) — a difference of two near-equal cumulative sums
                // that collapses to ~0 and flips sign in the quasi-isothermal upper
                // atmosphere, seeding a huge negative top radiation that back-substitution
                // smeared into NaN temperatures aloft (see project_multilayer_radiation).
                double aa, bb, cc, dd;
                for (int i = i_mount; i <= i_trop; i++) {
                    if (i == i_mount) {
                        aa = 0.0;
                        bb = -2.0 * m.radiation.x[i][j][k];
                        cc = m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k];
                        dd = -(1.0 - m.albedo.y[j][k]) * m.short_wave_radiation[j];
                        alfa[i] = -cc / bb;
                        beta[i] = +dd / bb;
                    }
                    if (i == i_mount + 1) {
                        aa = m.radiation.x[i - 1][j][k];
                        bb = -2.0 * m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];
                        cc = m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k];
                        dd = AA[i];
                        alfa[i] = -cc / (bb + aa * alfa[i - 1]);
                        beta[i] = +(dd - aa * beta[i - 1]) / (bb + aa * alfa[i - 1]);
                    }
                    if (i == i_mount + 2) {
                        aa = m.epsilon.x[i - 1][j][k] * m.radiation.x[i - 1][j][k];
                        bb = -2.0 * m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];
                        cc = m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k];
                        dd = -AA[i - 1] + AA[i] + CC[i - 1][i];
                        alfa[i] = -cc / (bb + aa * alfa[i - 1]);
                        beta[i] = +(dd - aa * beta[i - 1]) / (bb + aa * alfa[i - 1]);
                    }
                    if (i > i_mount + 2) {
                        aa = m.epsilon.x[i - 1][j][k] * m.radiation.x[i - 1][j][k];
                        bb = -2.0 * m.epsilon.x[i][j][k] * m.radiation.x[i][j][k];
                        // Top row (i == i_trop) has no layer above -> cc = 0. Guarding this
                        // also avoids the out-of-bounds read of radiation[i_trop+1] that the
                        // extended sweep would otherwise make.
                        cc = (i < i_trop) ? m.epsilon.x[i + 1][j][k] * m.radiation.x[i + 1][j][k]
                                          : 0.0;
                        dd = -AA[i - 1] + AA[i] - CA[i - 1] + CA[i];
                        alfa[i] = -cc / (bb + aa * alfa[i - 1]);
                        beta[i] = +(dd - aa * beta[i - 1]) / (bb + aa * alfa[i - 1]);   // FIX: aa (sub-diagonal), was alfa[i]
                    }
                }

                // Back-substitution (Thomas). The top unknown is beta[i_trop] because the
                // top row's alfa[i_trop] = 0 (cc = 0 — no layer above).
                m.radiation.x[i_trop][j][k] = beta[i_trop];
                for (int i = i_trop - 1; i >= 0; i--)
                    m.radiation.x[i][j][k] = alfa[i] * m.radiation.x[i + 1][j][k] + beta[i];

                // Radiation -> temperature (add back the reference emission, invert sigma T^4).
                for (int i = 0; i <= i_trop; i++) {
                    m.radiation.x[i][j][k] = radiation_original[i] + m.radiation.x[i][j][k];
                    m.t.x[i][j][k] = pow(m.radiation.x[i][j][k] / m.sigma, 0.25) / m.t_0;
                }

                // ---- Surface energy balance (radiative-CONVECTIVE) ----
                // The tridiagonal inversion left the surface T ~insensitive to the longwave
                // opacity, so CO2 did not warm the surface. Replace the surface value with an
                // explicit balance:  (1-albedo)*SW + L_down = sigma*T_s^4 + c_H*(T_s - T_air1).
                //  - L_down = downwelling longwave (back-radiation) = sum of atmospheric-layer
                //    emissions transmitted down to the surface; it RISES with CO2/H2O emissivity,
                //    so more greenhouse -> warmer surface (the physically-correct response).
                //  - c_H*(T_s - T_air1): bulk turbulent (sensible+latent) flux to the lowest air
                //    layer, which keeps the surface off the pure-radiative overheating —
                //    radiative-convective, not pure radiative. sigma*T_s^4 is linearised about
                //    the current surface T_s0 (one Newton step). See project_multilayer_radiation.
                double L_down = 0.0, trans = 1.0;
                for (int i = i_mount + 1; i <= i_trop; i++) {
                    L_down += m.epsilon.x[i][j][k] * m.sigma
                            * pow(m.t.x[i][j][k] * m.t_0, 4.0) * trans;     // layer i emission reaching surface
                    trans  *= (1.0 - m.epsilon.x[i][j][k]);                // attenuation through layer i
                }
                const double SW_abs = (1.0 - m.albedo.y[j][k]) * m.short_wave_radiation[j];
                const double T_air1 = m.t.x[i_mount + 1][j][k] * m.t_0;     // lowest air-layer T [K]
                const double T_s0   = m.t.x[i_mount][j][k] * m.t_0;         // linearisation point [K]
                const double c_H    = 15.0;                                // bulk turbulent transfer [W/m2/K]
                const double dsigT4 = 4.0 * m.sigma * T_s0 * T_s0 * T_s0;
                // ATHAD: the surface is molten, so it supplies heat from below as well as
                // absorbing it from above. A quenching magma ocean radiates far more than the
                // modern Earth's 0.09 W/m2, and at 1500 K this term is plausibly comparable to
                // the absorbed solar — omitting it would let the surface cool as if it were
                // rock. It enters the balance exactly as absorbed shortwave does.
                const double T_s    = (SW_abs + m.geothermal_flux + L_down
                                       - m.sigma * pow(T_s0, 4.0)
                                       + dsigT4 * T_s0 + c_H * T_air1) / (dsigT4 + c_H);
                m.radiation.x[i_mount][j][k] = m.sigma * pow(T_s, 4.0);
                m.t.x[i_mount][j][k]         = T_s / m.t_0;

                // De-kink the surface radiative step in the DIAGNOSTIC radiation profile only.
                // The 1-point surface energy balance (sigma T_s^4 at i_mount) and the column
                // Thomas solve above it are computed separately, leaving a sharp discontinuity
                // at the surface. One light 1-2-1 pass over the lowest layers softens it for
                // plotting. This touches ONLY radiation.x — t.x / the surface-balance T_s (and
                // hence the CO2 surface sensitivity) are left exactly as computed above.
                {
                    const int i_top_sm = std::min(i_mount + 4, i_trop);
                    double r_orig[5];                            // originals (i_mount .. i_top_sm)
                    for (int i = i_mount; i <= i_top_sm; i++)
                        r_orig[i - i_mount] = m.radiation.x[i][j][k];
                    for (int i = i_mount + 1; i < i_top_sm; i++)  // interior 1-2-1
                        m.radiation.x[i][j][k] = 0.25 * r_orig[i - 1 - i_mount]
                                               + 0.5  * r_orig[i - i_mount]
                                               + 0.25 * r_orig[i + 1 - i_mount];
                    if (i_top_sm > i_mount)                       // surface: one-sided blend toward air
                        m.radiation.x[i_mount][j][k] = 0.5 * r_orig[0] + 0.5 * r_orig[1];
                }
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for MultiLayerRadiation\n", elapsed.count() * 1e-9);
        cout << "      RadiationMultiLayer ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
