#pragma once

#include "cAtmosphereModel.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace AtomUtils;


class SaturationAdjustment {
public:
    explicit SaturationAdjustment(cAtmosphereModel& model)
        : m(model)
    {}


    void run() {
        using namespace std;

        cout << endl << endl << endl << "      SaturationAdjustment" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        cout.precision(9);

        computeSteps();
        adjustSaturation();
        applyTopography();
        clampAndFade();
        printReport();

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for SaturationAdjustment\n", elapsed.count() * 1e-9);

        cout << "      SaturationAdjustment ended" << endl;
    }

private:
    cAtmosphereModel& m;

    std::vector<double> step;

    // Precomputed constants
    static constexpr double fade_K        = 5.0;                        // transition half-width in Kelvin
    static constexpr int    iter_prec_end = 20;

    void computeSteps() {
        step.resize(m.im);
        for (int i = 0; i < m.im; i++)
            step[i] = m.get_layer_height(i+1) - m.get_layer_height(i);
    }

    void adjustSaturation() {
        const double inv_t_0      = 1.0 / m.t_0;
        const double t_range_inv  = 1.0 / (m.t_0 - m.t_00);
        const double lv_over_cp   = m.lv / m.cp_l;
        const double ls_over_cp   = m.ls / m.cp_l;

        // Surface row is skipped below, but its condensation source must still be
        // cleared every call: S_c_c.x[0] feeds the ice schemes' cloud-water source,
        // and over ocean (i_mount == 0) they copy S_c_c.x[0] onto itself, so a stale
        // nonzero value would silently re-condense surface cloud and defeat the skip.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++)
            for (int k = 0; k < m.km; k++)
                m.S_c_c.x[0][j][k] = 0.0;

        // Start at i = 1: the surface row (i = 0) is a vapour SOURCE only and must
        // never be saturation-adjusted in place. waterVapourEvaporation injects c[0]
        // from the warm ocean toward c_eq ~ q_sat every moist iter; if adjustSaturation
        // then condenses that just-injected vapour into cloud at the same cell, cloud
        // water doubles every ~3-4 iters (Cook Inlet runaway: cloud 5e-13 -> 32 over
        // iters 315-363). Once cloud + ice > 1 the r_humid denominator
        // (1 + 0.608*c - cloud - ice) flips negative, -grad(p)/rho inverts, NaN cascade.
        // Physically, condensation only happens once a parcel has risen, cooled
        // adiabatically and reached its LCL aloft; the warm ocean surface is below
        // saturation by definition. So cloud forms from i = 1 up, not at i = 0.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < m.im - 1; i++) {
            for (int j = 0; j < m.jm; j++) {
                double *S_c_c_row = m.S_c_c.x[i][j];
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row     = m.t.x[i][j];
                double *p_row     = m.p_stat.x[i][j];
                double  dt_dim    = step[i] / 1.6;

                for (int k = 0; k < m.km; k++) {
                    S_c_c_row[k] = 0.0;

                    double q_v_old = std::max(0.0, c_row[k]);
                    double q_c_old = std::max(0.0, cloud_row[k]);
                    double q_i_old = std::max(0.0, ice_row[k]);

                    double T       = t_row[k] * m.t_0;
                    // Cap T to the Magnus-formula validity range. Above ~101°C the
                    // saturation vapor pressure exceeds p_local and the q_sat fallback
                    // (ep * 1e-5) collapses by 5000×, triggering runaway condensation
                    // that releases lv/cp * q_v ≈ 87 K of latent heat per call and drives
                    // T further out of range. Cap at 60°C — well above any physical
                    // surface temperature — so a single corrupt cell cannot poison the run.
                    constexpr double T_max = 333.15;
                    if (T > T_max) {
                        T = T_max;
                        t_row[k] = T_max / m.t_0;
                    }
                    double p_local = p_row[k];

                    double E_sat  = m.hp * exp_func(T, 17.2694, 35.86);
                    double q_sat  = (p_local > E_sat)
                        ? m.ep * E_sat / (p_local - E_sat)
                        : m.ep * 1e-5;

                    const double alpha_entry = 1.0 / (1.0 + std::exp(-(T - m.t_00) / fade_K));

                    if ((q_v_old > q_sat && alpha_entry > 0.01) ||
                        (q_v_old < q_sat &&
                        (q_c_old > 1e-12 || q_i_old > 1e-12))) {

                        double q_v_b   = q_v_old;
                        double q_c_b   = q_c_old;
                        double q_i_b   = q_i_old;
                        double q_v_hyp = q_sat;
                        const double T_original = t_row[k] * m.t_0;

                        for (int iter = 1; iter <= iter_prec_end; iter++) {
                            double CND = std::max(0.0, std::min(1.0,
                                (T - m.t_00) * t_range_inv));
                            double DEP = 1.0 - CND;

                            double d_q_v = q_v_hyp - q_v_b;

                            if (d_q_v > 0) {
                                double max_evap = q_c_b + q_i_b;
                                if (d_q_v > max_evap) d_q_v = max_evap;
                            }

                            double d_cnd = d_q_v * CND;
                            double d_dep = d_q_v * DEP;

                            q_v_b += d_q_v;
                            q_c_b  = std::max(0.0, q_c_b - d_cnd);
                            q_i_b  = std::max(0.0, q_i_b - d_dep);

                            T -= lv_over_cp * d_cnd + ls_over_cp * d_dep;

                            double E_sat = m.hp * exp_func(T, 17.2694, 35.86);
                            double E_Ice = m.hp * exp_func(T, 21.8746, 7.66);
                            double q_sat = (p_local > E_sat)
                                ? m.ep * E_sat / (p_local - E_sat)
                                : m.ep * 1e-5;
                            double q_Ice = (p_local > E_Ice)
                                ? m.ep * E_Ice / (p_local - E_Ice)
                                : m.ep * 1e-5;

                            double q_sum = q_c_b + q_i_b;
                            double q_v_target = (q_sum > 1e-12)
                                ? (q_c_b * q_sat + q_i_b * q_Ice) / q_sum
                                : ((T >= m.t_0) ? q_sat : q_Ice);

                            // Adaptive Newton-damped update. The previous fixed 0.5 under-
                            // relaxation is UNSTABLE when the latent-heat gain
                            // G = (L/cp)*dq_sat/dT exceeds 3 — i.e. warm cells (T>~25C,
                            // q_sat>~20 g/kg) where dq_sat/dT is steep: the fixed-point map
                            // derivative 0.5*(1-G) then has magnitude >1, so the loop
                            // OSCILLATES and never converges, leaving RH stuck at 120-139%
                            // with cloud present (thermodynamically impossible) and feeding
                            // the coastal precip runaway. omega = 1/(1+G) drives the map
                            // derivative 1-omega*(1+G) to 0 (stable, ~Newton-optimal) at all T.
                            // dq_sat/dT from Clausius-Clapeyron: q_sat*L/(Rv*T^2).
                            // project_overprecip_saturation_injection.
                            const double Rv = 461.5;                    // [J/(kg K)] water-vapour gas constant
                            const double inv_RvT2 = 1.0 / (Rv * T * T);
                            const double Gain = CND * lv_over_cp * (q_sat * m.lv * inv_RvT2)
                                              + DEP * ls_over_cp * (q_Ice * m.ls * inv_RvT2);
                            const double omega = 1.0 / (1.0 + Gain);
                            q_v_hyp = q_v_b + omega * (q_v_target - q_v_b);

                            if (fabs(q_v_b / q_v_hyp - 1.0) <= 1.0e-6)
                                break;
                        }

                        q_c_b = std::max(0.0, q_c_b);
                        q_i_b = std::max(0.0, q_i_b);

                        // Cap T after the Newton loop. The q_v_hyp = 0.5*(q_v_target + q_v_b)
                        // damping is too weak when dq_sat/dT is steep (marginal saturation),
                        // so the iteration's amplitude grows. Within one call T can swing from
                        // a physical entry value into the Magnus-cliff regime (>101 °C at
                        // p = 1080 hPa), and the write-back below would persist that bad value.
                        // The entry-time cap is not enough because the runaway happens during
                        // the loop, not between calls.
                        if (T > T_max) T = T_max;

                        if (!std::isnan(T) && !std::isnan(q_v_b)) {
                            S_c_c_row[k] = alpha_entry * (q_c_b - q_c_old) / dt_dim;
                            c_row[k]     = q_v_old + alpha_entry * (q_v_b - q_v_old);
                            cloud_row[k] = q_c_old + alpha_entry * (q_c_b - q_c_old);
                            ice_row[k]   = q_i_old + alpha_entry * (q_i_b - q_i_old);
                            t_row[k]     = (T_original + alpha_entry * (T - T_original)) * inv_t_0;
                        }
 
                        if (T < m.t_00) {
                            cloud_row[k] = 0.0;
                            ice_row[k]   = 0.0;
                        }
                    }
                }
            }
        }
    }

    void applyTopography() {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                int i_mount = m.i_topography[j][k];
                m.c.x[0][j][k]     = m.c.x[i_mount][j][k];
                m.cloud.x[0][j][k] = m.cloud.x[i_mount][j][k];
                m.ice.x[0][j][k]   = m.ice.x[i_mount][j][k];
                if (AtomUtils::is_finite_safe(m.t.x[i_mount][j][k]))
                    m.t.x[0][j][k] = m.t.x[i_mount][j][k];
            }
        }
    }

    void clampAndFade() {
        const double inv_t_0    = 1.0 / m.t_0;
        const double lv_over_cp = m.lv / m.cp_l;
        const double ls_over_cp = m.ls / m.cp_l;
        // Defensive physical bounds. T_max mirrors the Magnus-validity cap used in
        // adjustSaturation; cloud_cap is ~50× the largest physical cloud/ice mixing
        // ratio (a few g/kg), so it never clips a real cloud — it only stops a runaway.
        constexpr double T_max     = 333.15;   // 60 °C
        constexpr double cloud_cap = 0.05;     // kg/kg condensate ceiling

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int j = 0; j < m.jm; j++) {
                double *c_row     = m.c.x[i][j];
                double *cloud_row = m.cloud.x[i][j];
                double *ice_row   = m.ice.x[i][j];
                double *t_row_nd  = m.t.x[i][j];
                double *p_row     = m.p_stat.x[i][j];

                for (int k = 0; k < m.km; k++) {
                    if (c_row[k]     < 0.0) c_row[k]     = 0.0;
                    if (cloud_row[k] < 0.0) cloud_row[k] = 0.0;
                    if (ice_row[k]   < 0.0) ice_row[k]   = 0.0;

                    double T_dim = t_row_nd[k] * m.t_0;
                    // Upper temperature bound first (no air parcel exceeds ~60 °C); keeps
                    // q_sat below finite and bounds the latent-heat release that follows.
                    if (T_dim > T_max) { T_dim = T_max; t_row_nd[k] = T_max * inv_t_0; }

                    // ---- Always-on supersaturation removal (ROOT FIX) ----
                    // adjustSaturation scales its condensation by alpha_entry, so in cold
                    // air (alpha ≪ 1) it leaves q_v ≫ q_sat — and cells colder than
                    // ~−60 °C (alpha ≤ 0.01) it skips entirely. With no upper bound on q_v
                    // anywhere, that residual supersaturation accumulated unbounded at the
                    // orographic saturation level in the cold fade zone (NZ Southern Alps,
                    // i≈12, T≈−40 °C → q_v→6, cloud→129, t→1e5 °C, then buoyancy drove a
                    // vertical-velocity runaway). Remove supersaturation here for EVERY
                    // cell, independent of the alpha fade, conserving water (excess → cloud
                    // above freezing, → ice below) and energy (latent heat → T). In a clean
                    // run the per-call excess is small (physical); the T_max recap below is
                    // the backstop if a transient ever drives a large excess.
                    const double p_local = p_row[k];
                    const double E_sat = m.hp * exp_func(T_dim, 17.2694, 35.86);
                    const double q_sat = (p_local > E_sat)
                        ? m.ep * E_sat / (p_local - E_sat)
                        : m.ep * 1e-5;
                    if (c_row[k] > q_sat) {
                        const double excess = c_row[k] - q_sat;
                        c_row[k] = q_sat;
                        if (T_dim >= m.t_00) {
                            cloud_row[k] += excess;
                            T_dim        += lv_over_cp * excess;
                        } else {
                            ice_row[k]   += excess;
                            T_dim        += ls_over_cp * excess;
                        }
                        if (T_dim > T_max) T_dim = T_max;   // backstop on the latent release
                        t_row_nd[k] = T_dim * inv_t_0;
                    }

                    // ---- Condensate upper bounds (CAP SAFETY NET) ----
                    if (cloud_row[k] > cloud_cap) cloud_row[k] = cloud_cap;
                    if (ice_row[k]   > cloud_cap) ice_row[k]   = cloud_cap;

                    // Existing cold fade: smoothly dry moisture toward 0 below t_00.
                    const double alpha = 1.0 / (1.0 + std::exp(-(T_dim - m.t_00) / fade_K));
                    c_row[k]     *= alpha;
                    cloud_row[k] *= alpha;
                    ice_row[k]   *= alpha;
                }
            }
        }
    }

    void printReport() const {
        // Diagnostic variables are not updated in the parallel loops;
        // kept here to preserve the original output contract.
        const bool   satadjust = false;
        const int    iter_prec = 0;
        const int    i_sat = 0, j_sat = 0, k_sat = 0;
        const double saturation = 0.0;

        if (!satadjust)
            std::cout << "      no saturation of water vapour in SaturationAdjustment found"
                      << std::endl;
        else
            std::cout << "      saturation of water vapour in SaturationAdjustment found"
                      << std::endl
                      << "      iter_prec = " << iter_prec << std::endl
                      << "      i_sat = "  << i_sat
                      << "   j_sat = "     << j_sat
                      << "   k_sat = "     << k_sat
                      << "   height_sat[m] = " << m.get_layer_height(i_sat)
                      << "   saturation[g/kg] = " << saturation * 1e3 << std::endl;

        if (iter_prec >= iter_prec_end)
            std::cout << std::endl
                      << "      no convergent solution found in SaturationAdjustment"
                      << std::endl
                      << "      iter_prec_end = " << iter_prec_end << std::endl
                      << "      iter_prec = "     << iter_prec     << std::endl
                      << "      results see above" << std::endl;
    }
};
