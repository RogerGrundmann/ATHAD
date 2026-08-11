#pragma once

#include "cAtmosphereModel.h"
#include "IceSchemeCommon.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <cstdio>

using namespace AtomUtils;


namespace TwoCatIce {
    // ice particle parameters
    constexpr double N_i_0    = 1.0e2;                                  // 1/m3
    constexpr double m_i_0    = 1.0e-12;                                // kg
    constexpr double m_i_max  = 1.0e-9;                                 // kg
    constexpr double m_s_0    = 3.0e-9;                                 // kg

    // microphysical rate coefficients
    constexpr double c_i_dep  = 1.3e-5;                                 // m3/(s*kg^(1/3))  [formula: c_i_dep * N_i[1/m3] * m_i^(1/3)[kg^(1/3)] * supersaturation -> 1/s]
    // Cloud->rain autoconversion rate. Raised from the COSMO 4.0e-4 to 1.0e-3 as the
    // rain-side partner of the reduced snow riming (c_rim_snow): reducing riming frees
    // cloud water that must reach the ground as RAIN rather than accumulate as cloud, so
    // the cloud->rain drain is sped up to match. Calibrated so global-mean precip stays
    // at the baseline (~5.9 mm/d) while the snow fraction drops from ~45% to ~16% — the
    // q_c_crit reservoir threshold is kept, so only excess cloud drains (no return to the
    // immediate-rain-out over-precipitation). project_snow_overproduction.
    constexpr double c_c_au   = 1.0e-3;                                 // 1/s (COSMO original 4.0e-4)
    constexpr double q_c_crit = 5.0e-4;                                 // [kg/kg] Kessler autoconversion threshold (~0.5 g/kg): cloud must accumulate before raining (project_overprecip_saturation_injection)
    constexpr double c_i_au   = 1.0e-3;                                 // 1/s COSMO
    constexpr double c_ac     = 0.24;                                   // m2/kg
    constexpr double c_rim    = 18.6;                                   // m2/kg
    // Reduced snow-side riming coefficient. Un-reduced riming (S_rim = c_rim*cloud*Snow)
    // dominated snow production ~10x every other source and monopolised the proportional
    // cloud-water limiter, routing ~90% of cloud water into SNOW and starving the
    // autoconversion->rain pathway — giving an unphysical ~45-58% snow fraction (Earth
    // ~5-10%). Reducing it (paired with the faster c_c_au autoconversion) lets the freed
    // cloud water convert to RAIN instead, dropping the snow fraction to ~16% with global
    // precip preserved. Warm-side shedding (S_shed) keeps the full c_rim.
    // project_snow_overproduction.
    constexpr double c_rim_snow = 18.6 / 5.0;                          // m2/kg (reduced snow riming)
    constexpr double c_agg    = 10.3;                                   // m2/kg
    constexpr double c_i_cri  = 0.24;                                   // m2/kg  [same structure as c_ac: c_i_cri * ice * Rain^(7/9)]
    constexpr double c_r_cri  = 3.2e-5;                                 // m2     [divided by m_i[kg] in formula -> effective m2/kg]
    constexpr double a_ev     = 1.0e-3;                                 // m2/kg
    constexpr double b_ev     = 5.9;                                    // m2*s/kg
    constexpr double c_s_dep  = 1.8e-2;                                 // m2/kg
    constexpr double b_s_dep  = 12.3;                                   // m2*s/kg
    constexpr double c_s_melt = 8.43e-5;                                // (m2*s)/(K*kg)
    constexpr double b_s_melt = 12.05;                                  // m2*s/kg
    constexpr double c_r_frz  = 3.75e-2;                                // m2/(K*kg)

    // temperature thresholds
    constexpr double t_nuc   = 267.15;                                  // K  (-6  °C)
    constexpr double t_d     = 248.15;                                  // K  (-25 °C)
    constexpr double t_hn    = 236.15;                                  // K  (-37 °C)
    constexpr double t_r_frz = 271.15;                                  // K  (-2  °C)

    // iteration / diagnostics
//    constexpr int iter_prec_end = 2;                                    // COSMO proposes 2 iterations
    constexpr int iter_prec_end = 3;                                    
//    constexpr int i_check       = 10;                                   // check level of P_rain
    constexpr int i_check       = 1;                                    // check level of P_rain

    // precomputed exponent fractions
    constexpr double exp_1_6  = 1.0 / 6.0;
    constexpr double exp_1_3  = 1.0 / 3.0;
    constexpr double exp_2_3  = 2.0 / 3.0;
    constexpr double exp_8_13 = 8.0 / 13.0;
    constexpr double exp_5_26 = 5.0 / 26.0;
    constexpr double exp_3_2  = 3.0 / 2.0;
}
/*
*
*/
// Two-Category-Ice-Scheme, COSMO-module from the German Weather Forecast,
// resulting the precipitation distribution formed of rain and snow.
class TwoCatIceScheme {
public:
    explicit TwoCatIceScheme(cAtmosphereModel& model)
        : m(model)
    {}

    void run() {
        using namespace std;
        using namespace TwoCatIce;

        cout << endl << endl << endl << "      TwoCategoryIceScheme" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        initArrays();
        computeColumns();
        applyBoundaryConditions();
        applyTopography();

        printReport();

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for TwoCategoryIceScheme\n", elapsed.count() * 1e-9);

        cout << "      TwoCategoryIceScheme ended" << endl << endl;
    }

private:
    cAtmosphereModel& m;

    // diagnostic output state (set during computeColumns)
    bool rain = false;
    bool snow = false;
    double P_sat = 0.0;
    double P_Snow = 0.0;
    int i_rain = 0, j_radain = 0, k_rain = 0;
    int i_snow = 0, j_snow = 0,   k_snow = 0;
    double height_rain = 0.0;
    double height_snow = 0.0;


    // ==================== INIT ====================
    void initArrays() {
        using namespace std;

        // Enforce non-negativity of the water species and warm-start the falling
        // precip flux from the previous iteration's accumulation each call. This
        // was previously the else-arm of an `if(m.Ma != 0)` split whose paleo arm
        // hard-zeroed the precip seed instead — but the member Ma is never assigned
        // (always 0, see project_atm_member_ma_latent_bug), so this arm always ran.
        // The two arms are equivalent for the precip field anyway (computeColumns
        // rebuilds P_rain/P_snow from a top-down integration with an iter_prec
        // fixed point each call, so the seed washes out); the only real difference
        // is the c/cloud/ice >= 0 clamp below, which we want for every run. The
        // dead paleo arm is dropped so a future Ma-member fix cannot silently
        // resurrect it and lose this clamp on paleo runs.
        #pragma omp parallel for collapse(2)
        for(int k = 0; k < m.km; k++){
            for(int j = 0; j < m.jm; j++){
                #pragma omp simd
                for(int i = 0; i < m.im; i++){
                    m.c.x[i][j][k]     = std::max(0.0, m.c.x[i][j][k]);
                    m.cloud.x[i][j][k] = std::max(0.0, m.cloud.x[i][j][k]);
                    m.ice.x[i][j][k]   = std::max(0.0, m.ice.x[i][j][k]);

                    m.P_rainn.x[i][j][k] = std::max(0.0, m.P_rainn.x[i][j][k]);
                    m.P_snown.x[i][j][k] = std::max(0.0, m.P_snown.x[i][j][k]);

                    m.P_rain.x[i][j][k] = m.P_rainn.x[i][j][k];
                    m.P_snow.x[i][j][k] = m.P_snown.x[i][j][k];
                }
            }
        }
    }
/*
*
*/
// ==================== COLUMN COMPUTATION ====================
    void computeColumns() {
        using namespace std;
        using namespace TwoCatIce;

        #pragma omp parallel for collapse(2)
        for(int k = 1; k < m.km-1; k++){
            for(int j = 1; j < m.jm-1; j++){

                double Rain_check = m.P_rain.x[i_check][j][k];

                // Column-private variables
                double q_sat = 0.0, E_sat = 0.0;
                double q_Ice = 0.0, E_Ice = 0.0;
                double dt_snow_dim = 0.0, dt_rain_dim = 0.0;
                double t_u, p_u;
                double m_i = m_i_max;

                double N_i = 0.0, S_nuc = 0.0, S_c_frz = 0.0, S_i_dep = 0.0,
                    S_c_au = 0.0, S_i_au = 0.0, S_d_au = 0.0,
                    S_ac = 0.0, S_rim = 0.0, S_shed = 0.0;
                double S_agg = 0.0, S_i_cri = 0.0, S_r_cri = 0.0, S_ev = 0.0,
                    S_s_dep = 0.0, S_i_melt = 0.0, S_s_melt = 0.0, S_r_frz = 0.0;

                // Work array for shell thickness (private per thread)
                double step[cAtmosphereModel::im];

                double P_rain_diff = 0.0;

                for(int iter_prec = 1; iter_prec <= iter_prec_end; iter_prec++){

                    m.P_rain.x[m.im-1][j][k] = 0.0;
                    m.P_snow.x[m.im-1][j][k] = 0.0;

                    for(int i = m.im-2; i >= 0; i--){

                        // Sub-terrain guard: cells with i < i_topography are inside the
                        // mountain. Their t/p/r_humid/cloud/ice values are sub-terrain copies
                        // (via subTerrainFill or similar paths) and produce unphysical S_c_au
                        // + S_ac + … → dP_rain that drove the iter-323 P_rain runaway at
                        // (i=6, j=30, k=209) — Gulf of Alaska / Cook Inlet inside-mountain
                        // cell. Zero rain/snow flux here so the downward integration only
                        // accumulates above the surface.
                        if (i < m.i_topography[j][k]) {
                            m.P_rain.x[i][j][k]        = 0.0;
                            m.P_snow.x[i][j][k]        = 0.0;
                            m.Precipitation.x[i][j][k] = 0.0;
                            m.S_v.x[i][j][k] = 0.0;
                            m.S_c.x[i][j][k] = 0.0;
                            m.S_i.x[i][j][k] = 0.0;
                            m.S_r.x[i][j][k] = 0.0;
                            m.S_s.x[i][j][k] = 0.0;
                            continue;
                        }

                        double Rain = m.P_rain.x[i][j][k];
                        double Snow = m.P_snow.x[i][j][k];

                        // Compute rain power terms only when raining
                        double Rain_pow_4_9  = 0.0;
                        double Rain_pow_7_9  = 0.0;
                        double Rain_pow_13_9 = 0.0;

                        if (Rain > 1e-12){
                            double rain_base  = pow(Rain, 1.0/9.0);     // ^1/9
                            double rain_pow_2 = rain_base * rain_base;  // ^2/9
                            double rain_pow_4 = rain_pow_2 * rain_pow_2;// ^4/9

                            Rain_pow_4_9  = rain_pow_4;
                            Rain_pow_7_9  = Rain/rain_pow_2;            // 1 - 2/9 = 7/9
                            Rain_pow_13_9 = Rain * rain_pow_4;          // 1 + 4/9 = 13/9
                        }

                        t_u = m.t.x[i][j][k] * m.t_0;                   // in K
                        p_u = m.p_stat.x[i][j][k];                      // in hPa

                        E_sat = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86); // saturation water vapour pressure for the water phase at t > 0°C in hPa
                        q_sat = m.ep * E_sat/(p_u - E_sat);             // relativ water vapour contents on ocean surface reduced by factor in kg/kg

                        E_Ice = m.hp * AtomUtils::exp_func(t_u, 21.8746, 7.66);
                        q_Ice = m.ep * E_Ice/(p_u - E_Ice);             // relativ water vapour contents on ocean surface reduced by factor in kg/kg

                        step[i] = m.get_layer_height(i+1)
                                - m.get_layer_height(i);                // local atmospheric shell thickness

                        double mass_layer = m.r_humid.x[i][j][k] * step[i]; // density * shell thickness

                        dt_rain_dim = step[i]/1.6;                      // adjusted rain fall time step by fixed velocities == 1.6 m/s
                        dt_snow_dim = step[i]/0.96;                     // adjusted snow fall time step by fixed velocities == 0.96 m/s


                        // number density of ice particles and mass of cloud ice
                        if((!(t_u > m.t_0)&&(!(t_u <= t_hn)))){
                            N_i = N_i_0 * exp(0.2 * (m.t_0 - t_u));     // in 1/m³

                            m_i = min((m.r_humid.x[i][j][k] * m.ice.x[i][j][k])
                                      /N_i, m_i_max);                   // in kg

                            if(m_i > m_i_max) m_i = m_i_max;            // m_i_max = 1.0e-9, in kg
                            if(m_i < m_i_0)   m_i = m_i_0;              // m_i_0 = 1.0e-12, in kg
                        }


                        // heterogeneous nucleation and depositional growth of cloud ice
                        if(m.ice.x[i][j][k] == 0.0){
                            if(((t_u < t_d)&&(m.c.x[i][j][k] >= q_Ice))
                                ||(((t_d <= t_u)&&(t_u <= t_nuc))
                                &&(m.c.x[i][j][k] >= q_sat)))

                                S_nuc = m_i_0/(m.r_humid.x[i][j][k]
                                        * dt_snow_dim) * N_i;           // nucleation of cloud ice, < I > in kg/(kg*s)
                        }else   S_nuc = 0.0;


                        // nucleation of cloud ice due to freezing of cloud water
                        if(t_u < t_hn)                                  // happens only below -37°C
                            S_c_frz = m.cloud.x[i][j][k]/dt_snow_dim;   // nucleation of cloud ice due to freezing of cloud water, < II > in kg/(kg*s)
                        else  S_c_frz = 0.0;


                        // deposition growth and sublimation of cloud ice
                        if(m.c.x[i][j][k] > q_Ice)                      // supersaturation
                            S_i_dep = c_i_dep * N_i * pow(m_i, exp_1_3) // c_i_dep = 1.3e-5, in m³/(s*kg)
                                * (m.c.x[i][j][k] - q_Ice);             // supersaturation, < III > in kg/(kg*s)

                        if(m.c.x[i][j][k] < q_Ice)                      // subsaturation, < III >
                            S_i_dep = max((- m.ice.x[i][j][k]/dt_snow_dim),
                                ((m.c.x[i][j][k] - q_Ice)/dt_snow_dim));// subsaturation, < III > in kg/(kg*s)


                        // autoconversion of cloud water to form rain (Kessler with threshold):
                        // cloud must exceed q_c_crit (~0.5 g/kg) before any rain forms. Without
                        // the threshold ANY cloud autoconverted instantly at c_c_au, so saturated
                        // columns rained out the moment they condensed (cloud~0 yet huge P_rain)
                        // -> gross precip ~30x NASA. See project_overprecip_saturation_injection.
                        if(m.cloud.x[i][j][k] > q_c_crit)               // c_c_au = 4.0e-4, in 1/s
                            S_c_au = c_c_au * (m.cloud.x[i][j][k] - q_c_crit); // cloud water to rain, cloud droplet collection, < IV > in kg/(kg*s)
                        else  S_c_au = 0.0;


                        // autoconversion of cloud ice to form snow due to aggregation
                        if(m.ice.x[i][j][k] > 0.0)                      // c_i_au = 1.0e-3, in 1/s
                            S_i_au = max(c_i_au * m.ice.x[i][j][k], 0.0); // cloud ice to snow, cloud ice crystal aggregation, < V > in kg/(kg*s)
                        else  S_i_au = 0.0;


                        // autoconversion of cloud ice to form snow due to deposition
                        if(S_i_dep > 0.0)
                            S_d_au = S_i_dep/(1.5 * (pow((m_s_0/m_i),
                                 exp_2_3) - 1.0));                      // autoconversion due to depositional growth of cloud ice, < VI > in kg/(kg*s)
                        else  S_d_au = 0.0;


                        // accretion of cloud water by raindrops (Kessler with the same
                        // q_c_crit threshold as autoconversion above): only the cloud water
                        // in excess of q_c_crit is collected into falling rain. Without the
                        // threshold accretion (c_ac=0.24, ~0.35x autoconv per unit cloud at
                        // 6 mm/d, growing with Rain) scavenged tenuous clouds the moment any
                        // rain advected/fell through them -> the broad ~2x over-precipitation
                        // that remained after the autoconversion fix. Keeping a q_c_crit cloud
                        // reservoir un-rainable makes accretion consistent with autoconversion.
                        // See project_overprecip_saturation_injection.
                        if((t_u >= m.t_0)&&(m.cloud.x[i][j][k] > q_c_crit))
                            S_ac = c_ac * (m.cloud.x[i][j][k] - q_c_crit) // c_ac = 0.24, in m²/kg, < VII > in kg/(kg*s)
                                   * Rain_pow_7_9;                      // in kg/(m² * s)
                        else  S_ac = 0.0;


                        // collection of cloud water by snow or graupel (riming)
                        if(t_u < m.t_0)
                            S_rim = c_rim_snow * m.cloud.x[i][j][k] * Snow;  // reduced riming (see c_rim_snow)
                        else  S_rim = 0.0;                              // riming rate of snow mass due to collection of supercooled cloud droplets, < VIII > in kg/(kg*s)

                        // collection of cloud water by wet snow to form rain (shedding)
                        if(t_u >= m.t_0)
                            S_shed = c_rim * m.cloud.x[i][j][k] * Snow; // c_rim = 18.6, m²/kg
                        else  S_shed = 0.0;                             // rate of water shed by melting wet snow particles, < IX > in kg/(kg*s)


                        // melting processes (snow/ice to rain at T > 0°C)
                        if(t_u > m.t_0){
                            // melting of falling snow
                            if(Snow > 1e-12){
                                S_s_melt = c_s_melt
                                    * (1.0 + b_s_melt * pow(Snow, exp_5_26))
                                    * (t_u - m.t_0) * pow(Snow, exp_1_3);
                                S_s_melt = std::min(S_s_melt, Snow/mass_layer); // Snow[kg/(m2*s)] / mass_layer[kg/m2] -> [1/s]
                            }else S_s_melt = 0.0;

                            // melting of cloud ice to cloud water/rain
                            S_i_melt = m.ice.x[i][j][k]/dt_snow_dim;
                        }else{
                            S_s_melt = 0.0;
                            S_i_melt = 0.0;
                        }


                        // collection of cloud ice by snow and rain (at T <= 0°C)
                        if(t_u <= m.t_0){
                            S_agg = c_agg * m.ice.x[i][j][k] * Snow;    // collection of cloud ice by snow particles, < X > in kg/(kg*s)

                            S_i_cri = c_i_cri * m.ice.x[i][j][k]        // c_i_cri = 0.24, m²/kg
                                * Rain_pow_7_9;                         // decrease in cloud ice due to collision/coalescence with raindrops, < XI > in kg/(kg*s)

                            S_r_cri = c_r_cri * m.ice.x[i][j][k]/m_i    // c_r_cri = 3.2e-5, m²/kg
                                      * Rain_pow_13_9;                  // decrease of rainwater due to freezing from ice crystal collection, < XII > in kg/(kg*s)
                        }else{
                            S_agg = 0.0;
                            S_i_cri = 0.0;
                            S_r_cri = 0.0;
                        }


                        // evaporation (rain to water vapour under subsaturation)
                        if(t_u > m.t_0 && Rain > 1e-12 && m.c.x[i][j][k] < q_sat){
                            S_ev = a_ev * (1.0 + b_ev * pow(Rain, exp_1_6))
                                   * (q_sat - m.c.x[i][j][k]) * Rain_pow_4_9;
                            S_ev = std::min(S_ev, Rain / mass_layer);   // Rain[kg/(m2*s)] / mass_layer[kg/m2] -> [1/s]
                        }else S_ev = 0.0;

                        // snow deposition/sublimation (at T < 0°C)
                        if(t_u < m.t_0 && Snow > 1e-12){
                            S_s_dep = c_s_dep * (1.0 + b_s_dep * pow(Snow, exp_5_26))
                                      * (m.c.x[i][j][k] - q_Ice) * pow(Snow, exp_8_13);
                            if(S_s_dep < 0.0)                           // sublimation limiting
                                S_s_dep = max(S_s_dep, -Snow/mass_layer); // Snow[kg/(m2*s)] / mass_layer[kg/m2] -> [1/s]
                        }else S_s_dep = 0.0;                            // no deposition above freezing or without snow

                        // freezing of rain to form snow
                        if(t_u <= m.t_0){                               // temperature below zero
                            double t_frz = max((t_r_frz - t_u), 0.0);

                            S_r_frz = c_r_frz * pow(t_frz, exp_3_2)     // immersion freezing and contact nucleation, < XVII > in kg/(kg*s)
                                * pow(Rain, exp_3_2);                   // c_r_frz = 3.75e-2, m²/(K*kg)
                        }else  S_r_frz = 0.0;


                        // cloud water limiter
                        double S_cloud_total = S_c_au + S_ac + S_rim + S_c_frz + S_shed;
                        double max_cloud_loss = m.cloud.x[i][j][k] / dt_snow_dim;

                        if(S_cloud_total > max_cloud_loss && S_cloud_total > 0){
                            double factor = max_cloud_loss/S_cloud_total;
                            S_c_au  *= factor;
                            S_ac    *= factor;
                            S_rim   *= factor;
                            S_c_frz *= factor;
                            S_shed  *= factor;
                        }

                        // ice consumption limiter: prevent over-depletion of cloud ice
                        double S_ice_total = S_i_au + S_d_au + S_agg + S_i_cri + S_i_melt;
                        double max_ice_loss = m.ice.x[i][j][k] / dt_snow_dim;

                        if(S_ice_total > max_ice_loss && S_ice_total > 0){
                            double factor = max_ice_loss / S_ice_total;
                            S_i_au  *= factor;
                            S_d_au  *= factor;
                            S_agg   *= factor;
                            S_i_cri *= factor;
                            S_i_melt *= factor;
                        }


                        // sinks and sources
                        m.S_v.x[i][j][k] = - m.S_c_c.x[i][j][k]         // sources and sinks of water vapour in kg/(kg*s)
                                           + S_ev - S_i_dep
                                           - S_s_dep - S_nuc;

                        m.S_c.x[i][j][k] =   m.S_c_c.x[i][j][k]         // sources and sinks of cloud water in kg/(kg*s)
                                           - S_c_au - S_ac
                                           - S_c_frz + S_i_melt
                                           - S_rim - S_shed;

                        m.S_i.x[i][j][k] =   S_nuc + S_c_frz            // sources and sinks of ice in kg/(kg*s)
                                           + S_i_dep - S_i_melt
                                           - S_i_au - S_d_au
                                           - S_agg - S_i_cri;

                        m.S_r.x[i][j][k] =   S_c_au + S_ac              // sources and sinks of rain in kg/(kg*s)
                                           - S_ev + S_shed
                                           - S_r_cri - S_r_frz
                                           + S_s_melt;

                        m.S_s.x[i][j][k] =   S_i_au + S_d_au            // sources and sinks of snow in kg/(kg*s)
                                           + S_agg + S_rim
                                           + S_s_dep + S_i_cri
                                           + S_r_cri + S_r_frz
                                           - S_s_melt;

                        // Per-cell precipitation cap (pragmatic bound on the gross-recirculation
                        // bursts). The stratiform condense<->recool recirculation produces a GROSS
                        // P_rain flux decoupled from the (few mm/d) net column water budget —
                        // Arabian/Scandinavia bursts ran to ~900-1024 mm/d and pulled the domain
                        // mean to ~10 mm/d (clean build). The old 0.1 mm/s (8640 mm/d) cap was ~3
                        // orders too loose to bind. Cap the rain & snow FLUX at every level to
                        // 50 mm/d so the accretion-amplified flux cannot accumulate beyond a
                        // physical ceiling; the total surface precip is additionally capped to
                        // 50 mm/d below. Band-aid (does not cure the recirculation) but bounds it
                        // deterministically. project_arabian_coast_precip_spike (option C).
                        constexpr double P_max = 50.0 / 86400.0;  // [kg/(m²s)] = [mm/s]  (= 50 mm/d)

                        // rain flux integration (top-down)
                        // S_i_melt excluded: cloud ice melts to cloud water (S_c), not
                        // directly to the falling rain flux; adding it here double-counts
                        // it with the S_c -> S_c_au -> dP_rain path.
                        double dP_rain = (S_c_au + S_ac + S_shed + S_s_melt
                            - S_ev - S_r_frz) * mass_layer;

                        m.P_rain.x[i][j][k] = std::min(P_max,
                            max(0.0, m.P_rain.x[i+1][j][k] + dP_rain));

                        // snow flux integration (top-down)
                        // S_i_melt excluded: suspended cloud ice and falling snow are
                        // separate categories; cloud ice melt does not remove mass from
                        // the falling snow flux.
                        double dP_snow = (S_i_au + S_d_au + S_rim + S_agg
                            + S_r_frz - S_s_melt) * mass_layer;

                        m.P_snow.x[i][j][k] = (t_u < m.t_0 && t_u >= m.t_000)
                            ? std::min(P_max, max(0.0, m.P_snow.x[i+1][j][k] + dP_snow))
                            : 0.0;

                        m.Precipitation.x[i][j][k] = std::min(P_max,
                            m.P_rain.x[i][j][k] + m.P_snow.x[i][j][k]);  // in mm/s, total capped at 50 mm/d

                    }  // end i


                    P_rain_diff = fabs(m.P_rain.x[i_check][j][k] - Rain_check) * 8.64e4;

                    if(P_rain_diff <= 1.0e-3){
                        std::cout.precision(10);
                        std::cout.setf(std::ios::fixed);
                        if((j == 90)&&(k == 180))  std::cout << endl
                            << " .... iteration ended .... no difference found" << endl
                            << " .... i_check = " << i_check
                            << " .... iter_prec = " << iter_prec << endl
                            << " .... temperature[C] = " << m.t.x[i_check][j][k] * m.t_0 - m.t_0
                            << " .... Rain_check = " << Rain_check * 8.64e4
                            << " .... P_rain[i_check] = " << m.P_rain.x[i_check][j][k] * 8.64e4
                            << " .... P_rain_diff = " << P_rain_diff << endl;
                        break;
                    }else{
                        std::cout.precision(10);
                        std::cout.setf(std::ios::fixed);
                        if((j == 90)&&(k == 180))  std::cout << endl
                            << " .... iteration ended .... differences are found" << endl
                            << " .... i_check = " << i_check
                            << " .... iter_prec = " << iter_prec << endl
                            << " .... temperature[C] = " << m.t.x[i_check][j][k] * m.t_0 - m.t_0
                            << " .... Rain_check = " << Rain_check * 8.64e4
                            << " .... P_rain[i_check] = " << m.P_rain.x[i_check][j][k] * 8.64e4
                            << " .... P_rain_diff = " << P_rain_diff << endl << endl << endl;
                    }
                }  // end iter_prec
            }  // end j
        }  // end k
    }
/*
*
*/
// ==================== BOUNDARY CONDITIONS ====================
    void applyBoundaryConditions() {
        IceSchemeCommon::extrapolateBC(m, m.P_rain, true);   // rain: phi-seam periodicity averaged
        IceSchemeCommon::extrapolateBC(m, m.P_snow, true);   // snow: phi-seam periodicity averaged
    }
/*
*
*/
// ==================== TOPOGRAPHY FILL ====================
    void applyTopography() {
        // Precipitation on mountains projected to sea level
        IceSchemeCommon::fillTopography(m, m.P_rain);
        IceSchemeCommon::fillTopography(m, m.P_snow);
    }
/*
*
*/
// ==================== DIAGNOSTIC REPORT ====================
    void printReport() const {
        using namespace std;

        if(!rain)
            cout << endl << "      no rain fall in TwoCategoryIceScheme found"
            << endl << endl;
        else
            cout << endl << "      rain fall in TwoCategoryIceScheme found"
            << endl
            << "      i_rain = " << i_rain
            << "      j_radain = " << j_radain
            << "   k_rain = " << k_rain
            << "   temperature[C] = " << m.t.x[i_rain][j_radain][k_rain] * m.t_0 - m.t_0
            << "   height_rain[m] = " << height_rain
            << "   P_sat[mm/d] = " << P_sat << endl << endl;

        if(!snow)
            cout << endl << "      no snow fall in TwoCategoryIceScheme found"
            << endl << endl;
        else
            cout << endl << "      snow fall in TwoCategoryIceScheme found"
            << endl
            << "      i_snow = " << i_snow
            << "      j_snow = " << j_snow
            << "   k_snow = " << k_snow
            << "   temperature[C] = " << m.t.x[i_snow][j_snow][k_snow] * m.t_0 - m.t_0
            << "   height_snow[m] = " << height_snow
            << "   P_Snow[mm/d] = " << P_Snow << endl << endl;
    }
};
/*
*
*/
