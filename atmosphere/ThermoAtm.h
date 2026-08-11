#pragma once

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

// Physical constants for precipitable water column integration
namespace PrecipWaterConstants {
    constexpr double HPA_TO_PA        = 100.0;
    constexpr double MIN_SAFE_TEMP    = 100.0;   // [K]
    constexpr double MIN_SAFE_PRESSURE = 1e-6;   // [hPa]
}

class ThermoAtm {
public:
    explicit ThermoAtm(cAtmosphereModel& model)
        : m(model)
    {}

    // ------------------------------------------------------------------
    void latentSensibleHeat()
    {
        using namespace std;
        cout << endl << endl << endl << "      LatentSensibleHeat" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double inv_2dr    = 1.0 / (2.0 * m.dr);
        const double inv_2dthe  = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi  = 1.0 / (2.0 * m.dphi);
        const double inv_sqrt3  = 1.0 / sqrt(3.0);

        std::vector<double> sinthe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            if (std::abs(sinthe_table[j]) < 1e-10)
                sinthe_table[j] = 1e-10;
        }

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                m.Q_Latent.x[0][j][k]   = 0.0;
                m.Q_Sensible.x[0][j][k] = 0.0;
            }
        }

        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 1; j < m.jm-1; j++) {
            for (int k = 1; k < m.km-1; k++) {

                double sinthe = sinthe_table[j];

                for (int i = 1; i < m.im-1; i++) {

                    if (is_land(m.h, i, j, k)) {
                        m.Q_Latent.x[i][j][k]   = 0.0;
                        m.Q_Sensible.x[i][j][k] = 0.0;
                        continue;
                    }

                    double rm           = m.rad.z[i];
                    double exp_rm       = 1.0 / (rm + 1.0);
                    double inv_rm       = 1.0 / rm;
                    double inv_rmsinthe = 1.0 / (rm * sinthe);

                    // Q_Latent is owned by RHS_Atm_Turb (the signed advective
                    // tendency × coeff_L written at every saturated cell). The
                    // earlier ThermoAtm path here wrote an unsigned magnitude
                    // diagnostic (lv · |∇q_sat|) every "moist" iter, racing the
                    // RHS writer on alternate iters and producing a 2Δt sawtooth
                    // in the field that confused diagnostics and (combined with
                    // the now-removed RHS_Atm_Turb.cpp:758 source term) drove the
                    // iter-358 NaN at 62°N upper-tropo.  See
                    // [[project-iter358-62n-dry-nan]] and
                    // [[project-precip-chain-fixes]].

                    // --- sensible heat ---
                    double dtdr   = (m.t.x[i+1][j][k] - m.t.x[i-1][j][k])
                                    * inv_2dr * exp_rm;
                    double dtdthe = (m.t.x[i][j+1][k] - m.t.x[i][j-1][k])
                                    * inv_2dthe * inv_rm;
                    double dtdphi = (m.t.x[i][j][k+1] - m.t.x[i][j][k-1])
                                    * inv_2dphi * inv_rmsinthe;

                    m.Q_Sensible.x[i][j][k] = m.cp_l * m.t_0 * inv_sqrt3
                        * sqrt(dtdr*dtdr + dtdthe*dtdthe + dtdphi*dtdphi);

                }  // i
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for LatentSensibleHeat\n", elapsed.count() * 1e-9);
        cout << "      LatentSensibleHeat ended" << endl;
    }

    // ------------------------------------------------------------------
    void waterVapourEvaporation()
    {
        using namespace std;
        cout << "\n\n\n      WaterVapourEvaporation  (model: " << m.evap_model << ")" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double conv_factor  = 8.64e4;                             // [mm/s] -> [mm/d]
        const double hPa_to_mmHg  = 0.750062;                           // [mmHg/hPa]
        const double ms_to_kmh    = 3.6;                                // [(km/h)/(m/s)]
        // Meyer (1915): E[mm/month] = C·(1 + W/16)·(e_s − e_a), with the deficit in
        // mmHg and C the open-water coefficient ≈ 11 (deep/large water, e.g. ocean) or
        // 15 (shallow ponds).  C is unit-identical for E in mm/month & deficit in mmHg
        // because inch→mm and inHg→mmHg share the 25.4 factor.  The previous value 0.36
        // was ~30× too small, giving ~0.18 mm/d (≈66 mm/yr) instead of the realistic
        // ~5 mm/d (≈1800 mm/yr) and far below the Dalton/Rohwer diagnostics.
        const double K_Meyer      = 11.0;                               // Meyer (1915) deep open-water coefficient [mm/month/mmHg]

        // Vertical spread: distribute c_eq over n_spread+1 levels with exp(-i) weights,
        // normalised so that sum_{i=0}^{n_spread} exp(-i) = 1 (moisture conserved).
//        const int    n_spread = 5;
        const int    n_spread = 3;
//        const int    n_spread = 2;
        const double r        = std::exp(-1.0);
        const double w_norm   = (1.0 - r) / (1.0 - std::pow(r, n_spread + 1)); // 1/sum

        // Determine which formula drives the surface-humidity (c) update.
        // All three are always computed for diagnostic output.
        enum class EvapModel { Dalton, Meyer, Rohwer } active;

        if      (m.evap_model == "Meyer")  active = EvapModel::Meyer;
        else if (m.evap_model == "Rohwer") active = EvapModel::Rohwer;
        else                               active = EvapModel::Dalton;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                if (is_land(m.h, 0, j, k)) {
                    m.Evaporation_Dalton.y[j][k] = 0.0;
                    m.Evaporation_Meyer.y[j][k]  = 0.0;
                    m.Evaporation_Rohwer.y[j][k] = 0.0;
                    m.Evaporation.y[j][k]        = 0.0;
                    continue;
                }

                double c_Dalton    = AtomUtils::C_Dalton(0, j, k, m.coeff_Dalton, m.u_0, m.u, m.v, m.w);
                                                                        // [mm/(h*hPa)]
                double p_stat_0jk  = m.p_stat.x[0][j][k];               // [hPa]
                // Floor the surface temperature at 180 K (below Earth's coldest-ever ~184 K, so
                // it never touches real physics) to backstop the marginal Gulf-of-Alaska crash:
                // an undamped i=0 coastal-surface T oscillation drove t.x[0] to ~160 K, making
                // E_sat/c_eq in the evaporation formula singular → c NaN. See
                // [[project_upper_velocity_secular_growth]].
                double t_u_base    = std::max(180.0, m.t.x[0][j][k] * m.t_0);  // [K], floored
                double precip_term = conv_factor * m.Precipitation.x[0][j][k]; // [mm/d]

                double vel_ms = sqrt((m.u.x[0][j][k] * m.u.x[0][j][k]
                                    + m.v.x[0][j][k] * m.v.x[0][j][k]
                                    + m.w.x[0][j][k] * m.w.x[0][j][k]) / 3.0) * m.u_0; // [m/s]
                double u_kmh  = vel_ms * ms_to_kmh;                     // [km/h]
                double p_mmHg = p_stat_0jk * hPa_to_mmHg;               // [mmHg]

                double E_sat  = (t_u_base >= m.t_0)                     // [hPa]
                    ? m.hp * AtomUtils::exp_func(t_u_base, 17.2694, 35.86)
                    : m.hp * AtomUtils::exp_func(t_u_base, 21.8746,  7.66);

                // Per-formula coefficients [mm/(d*hPa)]: E = coeff * sat_deficit_hPa
                //   Dalton: from wind-dependent mass-transfer coefficient
                //   Meyer (1915): K_M*(1+u/16)/30 [mm/month/mmHg -> mm/d/hPa]
                //   Rohwer (1931): 0.771*(1.465-0.000732*p)*(0.44+0.0733*u) [mm/d/mmHg -> mm/d/hPa]
                double coeff_D = c_Dalton * 24.0;
                double coeff_M = K_Meyer * hPa_to_mmHg * (1.0 + u_kmh / 16.0) / 30.0;
                double coeff_R = 0.771 * (1.465 - 0.000732 * p_mmHg)
                               * (0.44  + 0.0733   * u_kmh) * hPa_to_mmHg;

                // Calm conditions (zero wind): skip c update, Dalton = 0.
                // Meyer and Rohwer retain their still-air (u=0) terms.
                if (c_Dalton <= 0.0) {
                    double e  = m.c.x[0][j][k] * p_stat_0jk / m.ep;     // [hPa]  vapour pressure at c_Dalton = 0
                    double sd = std::max(0.0, E_sat - e);               // [hPa]  saturation deficit at c_Dalton = 0
                    m.Evaporation_Dalton.y[j][k] = 0.0;
                    m.Evaporation_Meyer.y[j][k]  = coeff_M * sd;        // [mm/d]
                    m.Evaporation_Rohwer.y[j][k] = coeff_R * sd;        // [mm/d]
                    m.Evaporation.y[j][k] = (active == EvapModel::Meyer)  ? m.Evaporation_Meyer.y[j][k]
                                          : (active == EvapModel::Rohwer) ? m.Evaporation_Rohwer.y[j][k]
                                          :                                 m.Evaporation_Dalton.y[j][k];
                    // Calm ocean cell: the wind-driven evaporation FLUX is zero (above), but the
                    // surface-layer vapour CONCENTRATION still equilibrates toward saturation (a
                    // windy no-precip cell converges to ~c_sat too). Relax c.x[0] toward c_sat with
                    // the same w_norm as the windy branch below, instead of leaving the
                    // RK4-unintegrated i=0 layer at a stale near-zero value — which showed as
                    // near-zero water-vapour patches over calm ocean (with normal "spots" at windy
                    // cells) in zonal cross-sections. Concentration only; the flux stays wind-limited.
                    double denom_calm = p_stat_0jk - (1.0 - m.ep) * E_sat;      // [hPa]
                    double c_sat_calm = (denom_calm > 0.0) ? m.ep * E_sat / denom_calm
                                                           : m.ep * E_sat / p_stat_0jk; // [kg/kg]
                    m.c_fix.y[j][k] = m.c.x[0][j][k];
                    m.c.x[0][j][k]  = m.c_fix.y[j][k] + (c_sat_calm - m.c_fix.y[j][k]) * w_norm;
                    continue;
                }

                // Saturation specific humidity (exact, not dilute approximation == Verdünnungsnäherung)
                double denom_guard = p_stat_0jk - (1.0 - m.ep) * E_sat; // [hPa]  guards that static pressure is greater than reduced current saturation pressure
                double c_sat       = (denom_guard > 0.0)
                    ? m.ep * E_sat / denom_guard                        // [kg/kg]
                    : m.ep * E_sat / p_stat_0jk;                        // dilute fallback

                // Active formula solves for c_eq balancing precipitation:
                //   E_active = coeff_active * (E_sat - c_eq*p/ep) = precip_term
                //   => c_eq = ep/p * (E_sat - precip_term / coeff_active)
                // If c_eq < 0, heavy rain saturates air -> clamp to c_sat.
                double coeff_active = (active == EvapModel::Meyer)  ?   coeff_M
                                    : (active == EvapModel::Rohwer) ?   coeff_R
                                    :                                   coeff_D;

                double c_eq = m.ep * (E_sat - precip_term / coeff_active) / p_stat_0jk;  // equilibrium humidity (dilute-approx saturation)
                c_eq = (c_eq < 0.0) ? c_sat : std::min(c_eq, c_sat);

                m.c_fix.y[j][k] = m.c.x[0][j][k];

                // Evaporations from current (pre-update) saturation deficit
                double e_cur  = m.c_fix.y[j][k] * p_stat_0jk / m.ep;    // [hPa]  current vapour pressure
                double sd_cur = std::max(0.0, E_sat - e_cur);           // [hPa]  current saturation deficit
 
                m.Evaporation_Dalton.y[j][k] = coeff_D * sd_cur;        // [mm/d]
                m.Evaporation_Meyer.y[j][k]  = coeff_M * sd_cur;        // [mm/d]
                m.Evaporation_Rohwer.y[j][k] = coeff_R * sd_cur;        // [mm/d]
                m.Evaporation.y[j][k] = (active == EvapModel::Meyer)  ? m.Evaporation_Meyer.y[j][k]
                                      : (active == EvapModel::Rohwer) ? m.Evaporation_Rohwer.y[j][k]
                                      :                                 m.Evaporation_Dalton.y[j][k];

                // c update: distribute c_eq across surface + n_spread levels above,
                // conserving total moisture (weights sum to 1).
                m.c.x[0][j][k] = m.c_fix.y[j][k] + (c_eq - m.c_fix.y[j][k]) * w_norm;  // i=0: weight = exp(0)*w_norm

                if (c_eq > 0.0) {
                    for (int i = 1; i <= n_spread; i++) {
                        double weight  = std::pow(r, i) * w_norm;
                        double t_i     = m.t.x[i][j][k] * m.t_0;
                        double p_i     = m.p_stat.x[i][j][k];
                        double E_i     = (t_i >= m.t_0)
                            ? m.hp * AtomUtils::exp_func(t_i, 17.2694, 35.86)
                            : m.hp * AtomUtils::exp_func(t_i, 21.8746,  7.66);
                        double denom_i = p_i - (1.0 - m.ep) * E_i;
                        double c_sat_i = (denom_i > 0.0) ? m.ep * E_i / denom_i : m.ep * E_i / p_i;
                        m.c.x[i][j][k] = std::min(m.c.x[i][j][k] + c_eq * weight, c_sat_i);
                    }
                }


/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if ((j == 90) && (k == 180)) cout << endl
                    << "  WaterVapourEvaporation" << endl
                    << "  j = " << j << "  k = " << k << endl
                    << "  p_stat_0jk = "     << p_stat_0jk
                    << "  p_mmHg = "         << p_mmHg << endl
                    << "  c_fix = "          << m.c_fix.y[j][k] * 1e3
                    << "  c_eq = "           << c_eq * 1e3
                    << "  c_sat = "          << c_sat * 1e3
                    << "  c = "              << m.c.x[0][j][k] * 1e3 << endl
                    << "  E_sat = "          << E_sat
                    << "  e_cur = "          << e_cur
                    << "  sd_cur = "         << sd_cur << endl
                    << "  c_Dalton = "       << c_Dalton
                    << "  coeff_D = "        << coeff_D
                    << "  coeff_M = "        << coeff_M
                    << "  coeff_R = "        << coeff_R << endl
                    << "  Evap_Dalton = "    << m.Evaporation_Dalton.y[j][k]
                    << "  Evap_Meyer = "     << m.Evaporation_Meyer.y[j][k]
                    << "  Evap_Rohwer = "    << m.Evaporation_Rohwer.y[j][k]
                    << "  Evap = "           << m.Evaporation.y[j][k]
                    << "  Prec = "           << precip_term << endl;
*/
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for WaterVapourEvaporation\n", elapsed.count() * 1e-9);
        cout << "      WaterVapourEvaporation ended" << endl;
    }

    // ------------------------------------------------------------------
    void standAtm_DewPoint_HumidRel()
    {
        using namespace std;
        cout << endl << endl << endl << "      StandAtm_DewPoint_HumidRel" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        const double lapse_rate = 6.5e-3;
        const double inv_ep     = 1.0 / m.ep;

        std::vector<int> i_trop_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            i_trop_table[j] = (j == 90)
                ? (int)m.tropopause_layers[j-1]
                : (int)m.tropopause_layers[j];
        }

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                int i_trop = i_trop_table[j];

                double TempStand_surface   = m.t.x[0][j][k] * m.t_0 - m.t_0;
                m.TempStand.x[0][j][k]    = TempStand_surface;

                for (int i = 0; i <= i_trop; i++) {

                    double t_u = m.t.x[i][j][k] * m.t_0;

                    double E = (t_u > m.t_0)
                        ? m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86)
                        : m.hp * AtomUtils::exp_func(t_u, 21.8746,  7.66);

                    double e           = m.c.x[i][j][k] * m.p_stat.x[i][j][k] * inv_ep;
                    bool   zero_vapour = (e == 0.0);
                    if (zero_vapour) e = 1e-3;

                    m.TempStand.x[i][j][k]    = TempStand_surface - lapse_rate * height_table[i];

                    double L    = std::log(e / m.hp);
                    double a_dp = (t_u > m.t_0) ? 21.8746 : 17.2694;
                    double b_dp = (t_u > m.t_0) ?  7.66   : 35.86;
                    m.TempDewPoint.x[i][j][k] = L * (m.t_0 - b_dp) / (a_dp - L);

                    m.HumidityRel.x[i][j][k]  = zero_vapour
                        ? 0.0
                        : std::min(e / E * 100.0, 100.0);

                }  // i

                // Above tropopause: copy tropopause values
                double ts_trop = m.TempStand.x[i_trop][j][k];
                double td_trop = m.TempDewPoint.x[i_trop][j][k];
                double hr_trop = m.HumidityRel.x[i_trop][j][k];

                for (int i = i_trop + 1; i < m.im; i++) {
                    m.TempStand.x[i][j][k]    = ts_trop;
                    m.TempDewPoint.x[i][j][k] = td_trop;
                    m.HumidityRel.x[i][j][k]  = hr_trop;
                }

                // Below-terrain fill: copy from topography surface downward
                int    i_mount   = m.i_topography[j][k];
                double ts_mount  = m.TempStand.x[i_mount][j][k];
                double td_mount  = m.TempDewPoint.x[i_mount][j][k];
                double hr_mount  = m.HumidityRel.x[i_mount][j][k];

                for (int i = i_mount - 1; i >= 0; i--) {
                    if (is_land(m.h, i, j, k)) {
                        m.TempStand.x[i][j][k]    = ts_mount;
                        m.TempDewPoint.x[i][j][k] = td_mount;
                        m.HumidityRel.x[i][j][k]  = hr_mount;
                    }
                }
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for StandAtm_DewPoint_HumidRel\n", elapsed.count() * 1e-9);
        cout << "      StandAtm_DewPoint_HumidRel ended" << endl;
    }

    // ------------------------------------------------------------------
    void forces()
    {
        using namespace std;
        cout << endl << endl << endl << "      ATOM: Forces" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double inv_2dr    = 1.0 / (2.0 * m.dr);
        const double inv_2dthe  = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi  = 1.0 / (2.0 * m.dphi);
        const double two_omega  = 2.0 * m.omega;
        const double omega2     = m.omega * m.omega;
        const double pres_coeff = 1e2 * m.p_0 / m.L_atm;
        const double r_Earth_m  = m.r_Earth * 1e3;

        std::vector<double> sinthe_table(m.jm);
        std::vector<double> costhe_table(m.jm);
        for (int j = 0; j < m.jm; j++) {
            sinthe_table[j] = sin(m.the.z[j]);
            costhe_table[j] = cos(m.the.z[j]);
        }

        #pragma omp parallel for collapse(2)
        for (int i = 1; i < m.im-1; i++) {
            for (int j = 1; j < m.jm-1; j++) {

                double rm           = m.rad.z[i];
                double exp_rm       = 1.0 / (rm + 1.0);
                double inv_rm       = 1.0 / rm;
                double sinthe       = sinthe_table[j];
                double costhe       = costhe_table[j];
                double inv_rmsinthe = 1.0 / (rm * sinthe);
                double abs_sinthe   = fabs(sinthe);

                double rad_dist  = (double)i * m.L_atm * exp_rm;
                double rad_Earth = rad_dist + r_Earth_m;

                for (int k = 1; k < m.km-1; k++) {

                    double w_ijk   = m.w.x[i][j][k] * m.u_0;
                    double Cor_r   = -two_omega * sinthe * w_ijk;
                    double Cor_the =  two_omega * costhe * w_ijk;
                    double Cor_phi =  two_omega * (-costhe * m.v.x[i][j][k]
                                                   + sinthe * m.u.x[i][j][k]) * m.u_0;

                    m.CoriolisForce.x[i][j][k] = m.Coriolis * m.r_air
                        * sqrt(Cor_r*Cor_r + Cor_the*Cor_the + Cor_phi*Cor_phi);

                    m.CentrifugalForce.x[i][j][k] =
                        m.centrifugal * m.r_air * omega2 * rad_Earth * (1.0 + abs_sinthe);

                    // Diagnostic buoyancy force for ParaView/Results — must match the
                    // perturbation-form body force applied in RHS_Atm.cpp (rhs_u), i.e.
                    // +g·ρ·(t − t̄(i)) referenced to the per-level mean t_ref_level[i].
                    const double t_ref_b = ((int)m.t_ref_level.size() == m.im)
                                           ? m.t_ref_level[i] : 1.0;
                    m.BuoyancyForce.x[i][j][k] = 1.0e-3 * m.buoyancy * m.r_humid.x[i][j][k] * m.g
                                                 * (m.t.x[i][j][k] - t_ref_b);

                    double dpdr   = (m.p_dyn.x[i+1][j][k] - m.p_dyn.x[i-1][j][k])
                                    * inv_2dr * exp_rm;
                    double dpdthe = (m.p_dyn.x[i][j+1][k] - m.p_dyn.x[i][j-1][k])
                                    * inv_2dthe * inv_rm;
                    double dpdphi = (m.p_dyn.x[i][j][k+1] - m.p_dyn.x[i][j][k-1])
                                    * inv_2dphi * inv_rmsinthe;

                    m.PresGradForce.x[i][j][k] =
                        - 1.0e-3 * sqrt(dpdr*dpdr + dpdthe*dpdthe + dpdphi*dpdphi) * pres_coeff;
                }  // k
            }  // j
        }  // i

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
        cout << "      ATOM: Forces ended" << endl;
    }

    // ------------------------------------------------------------------
    void adjustTemperatureIC(double** t, int jm, int km)
    {
        const double inv_t0 = 1.0 / m.t_0;
        const int    k_half = (km - 1) / 2;

        #pragma omp parallel for
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++)
                t[j][k] = (t[j][k] + m.t_0) * inv_t0;
            t[j][k_half] = (t[j][k_half + 1] + t[j][k_half - 1]) * 0.5;
            m.temperature_NASA.y[j][k_half] =
                (m.temperature_NASA.y[j][k_half + 1] +
                 m.temperature_NASA.y[j][k_half - 1]) * 0.5;
        }
    }

    // ------------------------------------------------------------------
    void precipitableWater()
    {
        using namespace std;
        using namespace PrecipWaterConstants;

        cout << "\n\n\n      PrecipitableWater" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.precipitable_water.y[j][k] = 0.0;

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                double column_sum = 0.0;

                for (int i = 0; i < m.im - 1; i++) {

                    const double t_actual  = m.t.x[i][j][k] * m.t_0;
                    const double p_actual  = m.p_stat.x[i][j][k];
                    const double q_mixing  = m.c.x[i][j][k];

                    if (t_actual < MIN_SAFE_TEMP || p_actual < MIN_SAFE_PRESSURE) {
                        m.PrecipitableWaterLocal.x[i][j][k] = 0.0;
                        continue;
                    }

                    const double e          = HPA_TO_PA * q_mixing * p_actual / m.ep;
                    const double a          = e / (m.R_WaterVapour * t_actual);
                    const double step       = m.get_layer_height(i + 1) - m.get_layer_height(i);
                    const double local_mass = a * step;

                    m.PrecipitableWaterLocal.x[i][j][k] = local_mass;
                    column_sum += local_mass / m.r_0_water;
                }

                m.precipitable_water.y[j][k] = column_sum * 1e3;   // m -> mm
            }  // k
        }  // j

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" Time measured: %.3f seconds for PrecipitableWater\n", elapsed.count() * 1e-9);
        cout << "      PrecipitableWater ended" << endl;
    }

    // ------------------------------------------------------------------
    void vegetationLand()
    {
        using namespace std;
        cout << "\n\n\n      VegetationLand" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        constexpr double max_vegetation_height_m = 4400.0;
        constexpr double min_vegetation_temp_c   = -40.0;

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        double max_Precipitation = 0.0;

        #pragma omp parallel for collapse(2) reduction(max:max_Precipitation)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                max_Precipitation = std::max(max_Precipitation, m.Precipitation.x[0][j][k]);

        double inv_max_Precipitation = 1.0 / std::max(max_Precipitation, 1e-5);

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];
                if (is_land(m.h, 0, j, k)
                    && height_table[i_mount] < max_vegetation_height_m
                    && m.t.x[i_mount][j][k] * m.t_0 - m.t_0 > min_vegetation_temp_c) {
                    m.Vegetation.y[j][k] = m.Precipitation.x[0][j][k] * inv_max_Precipitation;
                } else {
                    m.Vegetation.y[j][k] = 0.0;
                }
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        cout << " time measured: " << fixed << setprecision(3)
             << elapsed.count() * 1e-9 << " seconds for VegetationLand" << endl;
        cout << "      VegetationLand ended" << endl;
    }

    // ------------------------------------------------------------------
    void printDataAtm()
    {
        using namespace std;
        cout << endl << endl << endl << "      print_data_atm" << endl;
        cout << endl << " properties of the atmosphere at the surface: " << endl;

        double temperature_NASA_average, temperature_average, temperature_expected_average,
               precipitablewater_average, precipitation_average, precipitation_NASA_average,
               co2_average, Evaporation_average, Evaporation_Dalton_average,
               Evaporation_Meyer_average, Evaporation_Rohwer_average;

        #pragma omp parallel sections
        {
            #pragma omp section
                { temperature_NASA_average    = AtomUtils::GetMean_2D(m.jm, m.km, m.temperature_NASA); }
            #pragma omp section
                { temperature_average         = (AtomUtils::GetMean_3D(m.jm, m.km, m.t) - 1.0) * m.t_0; }
            #pragma omp section
                { precipitablewater_average   = AtomUtils::GetMean_2D(m.jm, m.km, m.precipitable_water); }
            #pragma omp section
                { precipitation_average       = 365.0 * AtomUtils::GetMean_3D(m.jm, m.km, m.Precipitation); }
            #pragma omp section
                { precipitation_NASA_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.precipitation_NASA); }
            #pragma omp section
                { co2_average                 = AtomUtils::GetMean_2D(m.jm, m.km, m.co2_total); }
            #pragma omp section
                { Evaporation_Dalton_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Dalton); }
            #pragma omp section
                { Evaporation_Meyer_average   = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Meyer); }
            #pragma omp section
                { Evaporation_Rohwer_average  = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation_Rohwer); }
            #pragma omp section
                { Evaporation_average         = 365.0 * AtomUtils::GetMean_2D(m.jm, m.km, m.Evaporation); }
        }
        temperature_expected_average =
            m.get_temperatures_from_curve(*m.get_current_time(), m.m_global_temperature_curve);

        cout.precision(2);
        cout << endl << endl;

        auto row = [](const char* n1, double v1, const char* u1,
                      const char* n2, double v2, const char* u2,
                      const char* n3, double v3, const char* u3) {
            cout << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n1 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v1 << setw(6) << u1 << "   "
                 << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n2 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v2 << setw(6) << u2 << "   "
                 << setiosflags(ios::left) << setw(40) << setfill('.')
                 << n3 << " = " << resetiosflags(ios::left) << setw(7)
                 << fixed << setfill(' ') << v3 << setw(6) << u3 << endl;
        };

        double precip_mm_a = precipitation_average * 8.64e4;

        row(" precipitable water average", precipitablewater_average, " mm",
            " precipitation average per year", precip_mm_a, " mm/a",
            " precipitation average per day", precip_mm_a / 365.0, " mm/d");

        // Per-component breakdown: which of P_rain / P_snow / P_graupel / P_conv
        // is climbing? Same units as the total: mean × 365·86400 gives mm/a.
        // Max scan reports the worst single cell with its (i,j,k) so a runaway
        // can be localised geographically.
        {
            const double s_per_year = 365.0 * 8.64e4;
            auto component = [&](const char* name, const Array& F){
                double mean_v = AtomUtils::GetMean_3D(m.jm, m.km, const_cast<Array&>(F));
                double mx = 0.0; int mi = 0, mj = 0, mk = 0;
                for(int i = 0; i < m.im; i++)
                    for(int j = 0; j < m.jm; j++)
                        for(int k = 0; k < m.km; k++){
                            double v = F.x[i][j][k];
                            if(v > mx){ mx = v; mi = i; mj = j; mk = k; }
                        }
                cout << "    " << setw(8) << setfill(' ') << name
                     << " mean = " << scientific << setprecision(3) << mean_v * s_per_year
                     << " mm/a   max = " << mx * s_per_year
                     << " mm/a  @(i=" << mi << ",j=" << mj
                     << ",k=" << mk << ")" << fixed << endl;
            };
            cout << " precipitation by component (mm/a, surface-equivalent):" << endl;
            component("Precip",    m.Precipitation);   // sanity check vs model total above
            component("P_rain",    m.P_rain);
            component("P_snow",    m.P_snow);
            component("P_graupel", m.P_graupel);
            component("P_conv",    m.P_conv);
        }

        row(" precipitable water average", precipitablewater_average, " mm",
            " precipitation NASA average per year", precipitation_NASA_average, " mm/a",
            " precipitation NASA average per day", precipitation_NASA_average / 365.0, " mm/d");

        row(" precipitable water average", precipitablewater_average, " mm",
            " Evaporation_Dalton_average per year", Evaporation_Dalton_average, " mm/a",
            " Evaporation_Dalton_average per day", Evaporation_Dalton_average / 365.0, " mm/d");

        row(" co2_average", co2_average, " ppm",
            " Evaporation_average per year", Evaporation_Dalton_average, " mm/a",
            " Evaporation_average per day", Evaporation_Dalton_average / 365.0, " mm/d");

        row(" co2_average", co2_average, " ppm",
            " Evaporation_Meyer_average per year", Evaporation_Meyer_average, " mm/a",
            " Evaporation_Meyer_average per day", Evaporation_Meyer_average / 365.0, " mm/d");

        row(" co2_average", co2_average, " ppm",
            " Evaporation_Rohwer_average per year", Evaporation_Rohwer_average, " mm/a",
            " Evaporation_Rohwer_average per day", Evaporation_Rohwer_average / 365.0, " mm/d");

        cout << endl;

        cout << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_NASA_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_NASA_average
             << setw(6) << " deg" << "   "
             << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_average
             << setw(6) << " deg" << "   "
             << setiosflags(ios::left) << setw(40) << setfill('.')
             << " temperature_expected_average" << " = " << resetiosflags(ios::left)
             << setw(7) << fixed << setfill(' ') << temperature_expected_average
             << setw(6) << " deg" << endl << endl << endl;

        cout << "      print_data_atm ended" << endl;
    }

    // ------------------------------------------------------------------
    void co2Atmosphere()
    {
        using namespace std;
        cout << endl << endl << endl << "      AGCM: co2_atmosphere" << endl;

        double t_paleo_add = 0.0;

        // The paleo CO2 increment is defined relative to the preceding slice. On a
        // single-Ma run there is no foregoing slice to difference against (and
        // get_previous_time() would throw), so leave t_paleo_add = 0 in that case.
        if ((!m.use_NASA_temperature || *m.get_current_time() > 0) && !m.is_first_time_slice())
            t_paleo_add =
                m.get_temperatures_from_curve(*m.get_current_time(),  m.m_global_temperature_curve)
              - m.get_temperatures_from_curve(*m.get_previous_time(), m.m_global_temperature_curve);

        m.co2_paleo = t_paleo_add * (3.2886 * t_paleo_add
            + 6.5772 * m.t_equat_modern - 32.8859);

        double co2_average = 3.2886 * m.t_equat_modern * m.t_equat_modern
            - 32.8859 * m.t_equat_modern + 102.2148;

        cout.precision(3);

        const char* temperature_comment = "      temperature increase at paleo times: ";
        const char* temperature_gain    = " t increase";
        const char* temperature_unit    = "°C ";
        const char* co_comment          = "      co2 increase at paleo times: ";
        const char* co_gain             = " co2 increase";
        const char* co_modern           = "      mean co2 at modern times: ";
        const char* co_paleo_str        = "      mean co2 at paleo times: ";
        const char* co_average_str      = " co2 modern";
        const char* co_average_pal      = " co2 paleo";
        const char* co_unit             = "ppm ";

        cout << setiosflags(ios::left) << setw(55) << setfill('.')
            << temperature_comment << resetiosflags(ios::left) << setw(13)
            << temperature_gain << " = " << setw(7) << setfill(' ')
            << t_paleo_add << setw(5) << temperature_unit << endl
            << setiosflags(ios::left) << setw(55) << setfill('.')
            << co_comment << resetiosflags(ios::left) << setw(12) << co_gain << " = "
            << setw(7) << setfill(' ') << m.co2_paleo << setw(5) << co_unit
            << endl << setw(55) << setfill('.') << setiosflags(ios::left) << co_modern
            << resetiosflags(ios::left) << setw(13) << co_average_str << " = "
            << setw(7) << setfill(' ') << co2_average << setw(5) << co_unit
            << endl << setw(55) << setfill('.') << setiosflags(ios::left)
            << co_paleo_str << resetiosflags(ios::left) << setw(13) << co_average_pal
            << " = " << setw(7) << setfill(' ') << co2_average + m.co2_paleo
            << setw(5) << co_unit << endl;

        double co2_max   = 397.0;
        double inv_h_top = 1.0 / m.get_layer_height(m.im - 1);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                double co2_surf      = m.t.x[0][j][k] * co2_max + m.co2_paleo;
                m.co2.x[0][j][k]    = co2_surf;

                int    i_mount       = m.i_topography[j][k];
                double co2_at_mount  = 0.0;

                if (i_mount > 0) {
                    double x_mount  = m.get_layer_height(i_mount) * inv_h_top;
                    co2_at_mount    = parabola_interp(m.co2_tropopause, co2_surf, x_mount);
                }

                for (int i = 1; i < m.im; i++) {
                    if (i <= i_mount && is_land(m.h, i, j, k)) {
                        m.co2.x[i][j][k] = co2_at_mount;
                    } else {
                        double x = m.get_layer_height(i) * inv_h_top;
                        m.co2.x[i][j][k] = parabola_interp(m.co2_tropopause, co2_surf, x);
                    }
                }

                if (i_mount > 0 && is_land(m.h, 0, j, k))
                    m.co2.x[0][j][k] = co2_at_mount;
            }
        }

        // CO2 sensitivity multiplier: scale the whole (ppm) field uniformly. co2_scale=1 leaves
        // the field as built; co2_scale=2 is a clean CO2 doubling (modern base + paleo increment),
        // so the mode-3 radiation sees 2x the ppm and the perturbation MLR(2x field) - MLR(280 ppm)
        // carries the doubling forcing. Leaves the paleo-CO2 construction (co2_max, co2_paleo)
        // untouched — this is a dedicated experiment knob (default 1.0).
        if (m.co2_scale != 1.0) {
            #pragma omp parallel for collapse(2) schedule(static)
            for (int j = 0; j < m.jm; j++)
                for (int k = 0; k < m.km; k++)
                    for (int i = 0; i < m.im; i++)
                        m.co2.x[i][j][k] *= m.co2_scale;
        }
        cout << "      AGCM: co2_atmosphere  co2_scale = " << m.co2_scale << endl;
        cout << "      AGCM: co2_atmosphere ended" << endl;
    }

    // ------------------------------------------------------------------
    void densities()
    {
        using namespace std;
        cout << "\n\n\n      PressureDensity" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        const double R_W_R_A_m1       = m.R_WaterVapour / m.R_Air - 1.0;
        const double inv_R_Air        = 1.0 / m.R_Air;
        const double scale            = 1e2 * inv_R_Air;
        const double p_sl_factor      = 1e-2 * m.r_air * m.R_Air;
        const double beta             = 42.0;                           // K, COSMO
        const double inv_beta         = 1.0 / beta;
        const double two_beta_g_inv_R = 2.0 * beta * m.g * inv_R_Air;

        std::vector<double> height_table(m.im);
        for (int i = 0; i < m.im; i++)
            height_table[i] = m.get_layer_height(i);

        #pragma omp parallel for collapse(2) schedule(dynamic, 4)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {

                // Floor surface T at 180 K (below Earth's coldest-ever ~184 K): coeff ∝ 1/T²,
                // so an unphysical cold surface (the coastal/Pamir T oscillation) blows coeff up
                // and drives the barometric sqrt negative. Backstops the sqrt guard at the
                // source. See [[project_upper_velocity_secular_growth]].
                const double t_u_0        = std::max(180.0, m.t.x[0][j][k] * m.t_0);
                const double p_sl         = p_sl_factor * t_u_0;
                const double t_0_inv_beta = t_u_0 * inv_beta;
                const double coeff        = two_beta_g_inv_R / (t_u_0 * t_u_0);

                for (int i = 0; i < m.im; i++) {
                    const double h_i = height_table[i];
                    const double t_u = m.t.x[i][j][k] * m.t_0;
                    // Guard the COSMO barometric sqrt: coeff = 2βg/(R_Air·t_u_0²) ∝ 1/T², so an
                    // anomalously cold surface column drives (1 - coeff·h_i) negative aloft →
                    // sqrt(NaN) → p_i/r_humid NaN → whole-field blow-up. The init-time copy of
                    // this formula (InitValues_Atm.cpp:644) already clamps with max(0,…); this
                    // in-loop version had dropped it. See [[project_upper_velocity_secular_growth]].
                    const double p_i = p_sl * exp(-t_0_inv_beta
                        * (1.0 - sqrt(std::max(0.0, 1.0 - coeff * h_i))));  // COSMO barometric formula

                    m.p_stat.x[i][j][k]  = p_i;
                    m.r_dry.x[i][j][k]   = scale * p_i / t_u;
                    m.r_humid.x[i][j][k] = scale * p_i
                        / ((1.0 + R_W_R_A_m1 * m.c.x[i][j][k]
                            - m.cloud.x[i][j][k] - m.ice.x[i][j][k]) * t_u);
                }

                const int    i_m = m.i_topography[j][k];
                m.p_stat_landscape.y[j][k] = m.p_stat.x[i_m][j][k];

                if (i_m >= 0 && i_m < m.im) {
                    m.p_stat.x[0][j][k]  = m.p_stat.x[i_m][j][k];
                    m.r_dry.x[0][j][k]   = m.r_dry.x[i_m][j][k];
                    m.r_humid.x[0][j][k] = m.r_humid.x[i_m][j][k];
                }
            }
        }

        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for PressureDensity\n", elapsed.count() * 1e-9);
        cout << "      PressureDensity ended" << endl;
    }

private:
    cAtmosphereModel& m;
};
