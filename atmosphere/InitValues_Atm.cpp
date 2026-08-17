#include "MixtureAtm.h"
#include "SaturationH2O.h"
#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;


// ============================================================================
// Water Vapor and Cloud Initialization Module - Final Improved Version
// ============================================================================

#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>

// ============================================================================
// Physical and Numerical Constants
// ============================================================================
namespace VaporCloudConstants {
    // Surface evaporation coefficients
    constexpr double COEFF_LAND = 0.74;                                 // Land surface evaporation coefficient
    constexpr double COEFF_OCEAN = 0.98;                                // Ocean surface evaporation coefficient
    
    // Moisture limits
    constexpr double Q_LIMIT = 3.0e-3;                                  // Maximum specific humidity [kg/kg]
    constexpr double Q_SCALING = 0.84;                                  // Empirical moisture scaling factor
    constexpr double FALLBACK_Q_FACTOR = 1e-5;                          // Fallback saturation factor

    // Relative humidity and cloud thresholds
    constexpr double RH_THRESHOLD = 85.0;                               // Cloud formation RH threshold [%]
    constexpr double CLOUD_SCALING = 0.0005;                            // Cloud density scaling factor
    
    // Stability parameters
    constexpr double LAPSE_RATE_REF = -0.0065;                          // Reference lapse rate [K/m]
    constexpr double STABILITY_SCALING = 500.0;                         // Stability weight scaling
    constexpr double STABILITY_IMPACT = 0.3;                            // Max stability influence
    
    // Dewpoint spread threshold
    constexpr double SPREAD_THRESHOLD = 2.0;                            // Temperature-dewpoint spread [K]

    // Magnus-Tetens formula coefficients
    // The Magnus coefficients that lived here are gone: every saturation calculation now
    // goes through SaturationH2O.h (IAPWS), which is valid across ATHAD's whole 1500 K
    // column rather than the ~320 K the Magnus fit was calibrated on.

    // Dewpoint calculation (inverse Magnus)
    constexpr double MAG_A = 17.27;
    constexpr double MAG_B = 237.3;
    constexpr double E0_HPA = 6.1078;                                   // Reference vapor pressure [hPa]
    
    // Safety limits
    constexpr double MIN_PRESSURE = 1e-10;
    constexpr double MIN_Q_SATUR = 1e-12;
}

using namespace VaporCloudConstants;

// ============================================================================
// Main Water Vapor and Cloud Initialization
// ============================================================================
void cAtmosphereModel::init_vapour_cloud() {                            // calculates initial water vapour and cloud distribution, no ice clouds are prepared
    std::cout << "\n\n\n      AGCM: init_vapour_cloud" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // Thread-local variables
    double p_u = 0.0;
    double t_u = 0.0;

    // Precompute inverse for efficiency
    const double inv_rh_range = 1.0 / (100.0 - RH_THRESHOLD);

    // ========================================================================
    // Main Computation Loop: Calculate Vapor and Cloud Fields
    // ========================================================================
    #pragma omp parallel for collapse(2) private(p_u, t_u)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {

            // Temperature thresholds (local constants for each thread)
            const double T_ice_end = t_000;
            const double T_freeze = t_0;

            const int i_mount = i_topography[j][k];

            // Process vertical column
            for (int i = 0; i < im-1; i++) {

                // Local variables for this grid point
                double E_sat_loc = 0.0;
                double q_sat_loc = 0.0;

                // Get local temperature and pressure
                t_u = t.x[i][j][k] * t_0;                               // [K]
                p_u = p_stat.x[i][j][k];                                // [hPa]

                // ------------------------------------------------------------
                // Calculate Saturation Vapor Pressure (Magnus-Tetens)
                // ------------------------------------------------------------
                const double E_wat = SaturationH2O::saturationPressure(t_u);
                const double E_ice = SaturationH2O::sublimationPressure(t_u);
                double E_sat;
                                                                                                                                                                                                        
                // Temperature-dependent phase transition
                if (t_u >= T_freeze) {
                    E_sat = E_wat;                                      // Pure liquid water
                } else if (t_u <= T_ice_end) {
                    E_sat = E_ice;                                      // Pure ice
                } else {
                    // Linear interpolation in mixed phase region
                    const double w = (t_u - T_ice_end) / (T_freeze - T_ice_end);
                    E_sat = w * E_wat + (1.0 - w) * E_ice;
                }                                                                                                                                                                                                     

                // Saturation MASS FRACTION, exact conversion. The dilute
                // ep*E/(p - E) with a fallback constant when E >= p is wrong in both
                // branches once water is the bulk gas — see SaturationH2O.h. This
                // routine is not called (see cAtmosphereModel::RunTimeSlice), but the
                // forbidden form left sitting in it is one uncomment away from being live.
                const double q_Satur = SaturationH2O::saturationMassFraction(
                    E_sat, p_u, AtmMixture::M_nonwater(c.x[i][j][k], co2.x[i][j][k],
                                                       m_comp.M_bg));

                // ------------------------------------------------------------
                // Surface Evaporation (only above topography)
                // ------------------------------------------------------------
                if (i >= i_mount) {
                    const double current_coeff = is_land(h, i, j, k) ? COEFF_LAND : COEFF_OCEAN;
                    E_sat_loc = E_sat * current_coeff;
                } else {
                    E_sat_loc = 0.0;
                }

                // Calculate specific humidity from evaporation
                if (E_sat_loc > 0.0) {
                    q_sat_loc = SaturationH2O::saturationMassFraction(
                        E_sat_loc, p_u, AtmMixture::M_nonwater(c.x[i][j][k],
                                                               co2.x[i][j][k], m_comp.M_bg));
                } else {
                    q_sat_loc = 0.0;
                }

                // ------------------------------------------------------------
                // Dewpoint Temperature and Spread
                // ------------------------------------------------------------
                const double e_actual   = (p_u * q_sat_loc) / (ep + q_sat_loc);
                const double log_val    = std::log(std::max(e_actual, MIN_PRESSURE) / E0_HPA);
//                const double t_dewpoint = (MAG_B * log_val) / (MAG_A - log_val) + 273.15;
//                const double spread     = t_u - t_dewpoint;

                // ------------------------------------------------------------
                // Atmospheric Stability Assessment
                // ------------------------------------------------------------
                // Bounds check: ensure i+1 is valid
                const double dT = (i < im-2) ? (t.x[i+1][j][k] - t.x[i][j][k]) : 0.0;
//                const double diff = dT - LAPSE_RATE_REF;

                // Stability weight: reduces moisture in stable conditions
//                const double stability_weight = 1.0 - STABILITY_IMPACT * std::tanh(diff * STABILITY_SCALING);

                // ------------------------------------------------------------
                // Cloud Formation Signal (based on relative humidity)
                // ------------------------------------------------------------
                const double rel_hum = (q_Satur > MIN_Q_SATUR) ? 
                                       (q_sat_loc / q_Satur * 100.0) : 0.0;
                double cloud_signal = (rel_hum - RH_THRESHOLD) * inv_rh_range;
                cloud_signal = std::clamp(cloud_signal, 0.0, 1.0);

                // ------------------------------------------------------------
                // Final Moisture Content (with stability and spread constraints)
                // ------------------------------------------------------------
//                const double q_final = Q_SCALING * q_sat_loc;

                // Apply moisture only when conditions are favorable:
                // - Small dewpoint spread (near saturation)
                // - Above topography
/*
                if (spread < SPREAD_THRESHOLD * stability_weight && i >= i_mount) {
                    const double q_weighted = q_final * stability_weight * cloud_signal;
                    c.x[i][j][k] = std::min(Q_LIMIT, q_weighted);
                } else {
                    c.x[i][j][k] = 0.0;
                }
*/

                const double RH_init = is_land(h, i, j, k) ? 0.60 : 0.75;
                c.x[i][j][k] = (i >= i_mount) ? RH_init * q_Satur : 0.0;                                                                                                                                         

                // Cloud density field
                cloud.x[i][j][k] = cloud_signal * CLOUD_SCALING;
                if (t_u < t_00)  cloud.x[i][j][k] = 0.0;
            }  // end i loop

            // Boundary conditions at top of domain
            c.x[im-1][j][k] = 0.0;
            cloud.x[im-1][j][k] = 0.0;

        }  // end k loop
    }  // end j loop

    // ========================================================================
    // Apply Surface Boundary Condition (copy topography values to surface)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
            
            // Bounds check before accessing array
            if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k)) {
                c.x[0][j][k] = c.x[i_mount][j][k];
            }
        }
    }
/*
    // ========================================================================                                                                                                             
    // Post-processing: smooth c and cloud horizontally
    // ========================================================================                                                                                                             
    {           
        // Precompute Gaussian weights for the stencil (constant for all levels/cells)
        const int R = 2;
        const int D = 2 * R + 1;

        std::vector<double> gauss_w(D * D);

        for (int dj = -R; dj <= R; dj++)
            for (int dk = -R; dk <= R; dk++)
                gauss_w[(dj + R) * D + (dk + R)] =
                    std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));


        #pragma omp parallel
        {
            std::vector<double> tmp(jm * km);       // thread-local scratch buffer

            auto smoothLevel = [&](auto& field, int i) {
                for (int j = 0; j < jm; j++)
                    for (int k = 0; k < km; k++)
                        tmp[j * km + k] = field.x[i][j][k];

                for (int j = 0; j < jm; j++) {
                    for (int k = 0; k < km; k++) {
                        if (i < i_topography[j][k]) continue;

                        double sum = 0.0;
                        double cnt = 0.0;

                        for (int dj = -R; dj <= R; dj++) {
                            const int jj = std::clamp(j + dj, 0, jm - 1);

                            for (int dk = -R; dk <= R; dk++) {
                                const int kk = (k + dk + km) % km;
//                                const double w = std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));   // if not use Gaussian weights
                                const double w = gauss_w[(dj + R) * D + (dk + R)];

                                if (i >= i_topography[jj][kk]) {
                                    sum += w * tmp[jj * km + kk];
                                    cnt += w;
                                }
                            }
                        }
                        field.x[i][j][k] = cnt > 0.0 ? sum / cnt : 0.0;
                    }
                }
            };

            // Each i-level is independent: threads take different levels
            #pragma omp for schedule(dynamic, 4)
            for (int i = 0; i < im - 1; i++) {     // im-1 stays 0 (top BC)
                smoothLevel(c,     i);
                smoothLevel(cloud, i);
            }

        } // end omp parallel

        // Re-apply surface BC
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const int i_mount = i_topography[j][k];

                if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k))
                    c.x[0][j][k] = c.x[i_mount][j][k];
            }
        }
    }
*/

    // ========================================================================
    // Performance Timing and Output
    // ========================================================================
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for init_vapour_cloud\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: init_vapour_cloud ended" << std::endl;
}
/*
*
*/
// ============================================================================
// Optional: Debug Output Function (compile with -DDEBUG_VAPOR_CLOUD)
// ============================================================================
#ifdef DEBUG_VAPOR_CLOUD

void cAtmosphereModel::debug_vapor_output(int i, int j, int k, 
                                          double t_u, double p_u, 
                                          double E_sat, double q_sat) {
    if ((j == 60) && (k == 87)) {
        std::cout.precision(5);
        std::cout.setf(std::ios::fixed);
        std::cout << "\n  WaterVapour Debug Output"
                  << "\n  Position: i=" << i << " j=" << j << " k=" << k
                  << "\n  p_static = " << p_u
                  << "\n  T = " << t_u << " K  (" << t_u - t_0 << " °C)"
                  << "\n  E_sat = " << E_sat
                  << "\n  0.84 * q_sat = " << q_sat
                  << "\n  Moisture c = " << c.x[i][j][k]
                  << "\n  Diff_c = " << c.x[i][j][k] - q_sat
                  << "\n" << std::endl;
    }
}
#endif
/*
*
*/
void cAtmosphereModel::init_tropopause_layers(){                                                                                                                                                         
    cout << endl << endl << endl << "      AGCM: init_tropopause_layers" << endl;                                                                                                                        
                                                                                                                                                                                                           

    int j_max = jm - 1;                                                                                                                                                                                  
    int j_half = j_max / 2;                                                                                                                                                                              
                  
    // Derive x_max so that Agnesi(tropopause_equator, x_max) == tropopause_pole exactly.                                                                                                                
    // Agnesi: a^3/(a^2+x^2) = b  =>  x = a * sqrt(a/b - 1)
    // Requires tropopause_equator > tropopause_pole (always true physically).                                                                                                                           
    double x_max = tropopause_equator                                                                                                                                                                    
                   * std::sqrt(tropopause_equator / tropopause_pole - 1.0);                                                                                                                              
                  
    // ATHAD: height -> level index must INVERT the stretched grid.
    //
    // What was here: round(h / L_atm). That is a level index only if the layers were
    // uniformly L_atm thick, and they are not — init_layer_heights builds
    //     height(i) = (exp(zeta * i / (im-1)) - 1) * L_atm
    // so L_atm is the AMPLITUDE of an exponential stretch, not a layer spacing (the same
    // confusion the metric-length comment in cAtmosphereModel.h warns about). The exact
    // inverse is
    //     i = (im-1) * ln(1 + h / L_atm) / zeta.
    //
    // The error is not small and it is not latent here. At the 300 km shell the pole's
    // 195 km convective top sits at level 52; the old formula returned round(195000/15719)
    // = 12, which is 13 km. VelocityInitializer builds the whole jet profile between the
    // surface and this level and applies a linear taper to zero from it to the domain top,
    // so the initial wind structure was compressed into the bottom 4 % of the atmosphere
    // and the remaining 96 % got the taper. On Earth (L_atm = 400 m, h = 11 km) it returned
    // 28 of 41 levels, which is wrong too — 4.9 km, not 11 — but wrong in a way that still
    // lands inside the troposphere, so it never showed.
    const double idx_scale = (double)(im - 1) / zeta;
    auto height_to_level = [&](double h) -> double {
        if(!(h > 0.0) || !(L_atm > 0.0) || !(zeta > 0.0)) return 0.0;
        const double i = round(idx_scale * std::log(1.0 + h / L_atm));
        return std::min(std::max(i, 0.0), (double)(im - 1));
    };
    // The forward map, written out rather than taken from get_layer_height(): this
    // routine runs in an omp section CONCURRENT with init_layer_heights(), so reading
    // m_layer_heights here would race the vector that fills it.
    auto level_to_height = [&](double i) {
        return (std::exp(zeta * i / (double)(im - 1)) - 1.0) * L_atm;
    };

    cout << "tropopause_pole=" << tropopause_pole
         << " x_max=" << x_max
         << " pole_index=" << height_to_level(tropopause_pole)
         << " (height " << level_to_height(height_to_level(tropopause_pole))
         << " m of " << level_to_height(im-1) << " m)" << endl;

    // Build symmetric cache of heights [m] and grid indices in one pass.
    std::vector<double> tropo_height_cache(jm);
    tropopause_layers = std::vector<double>(jm);

    for(int j = 0; j <= j_half; j++){
        double x = x_max * (double)(j_half - j) / (double)j_half;
        double h = AtomUtils::Agnesi(tropopause_equator, x);
        tropo_height_cache[j]       = h;
        tropo_height_cache[j_max-j] = h;
        tropopause_layers[j]        = height_to_level(h);
        tropopause_layers[j_max-j]  = tropopause_layers[j];
    }
                                                                                                                                                                                                           
    #pragma omp parallel for schedule(static)                                                                                                                                                            
   for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){                                                                                                                                                                     
            Tropopause.y[j][k] = tropo_height_cache[j];
        }                                                                                                                                                                                                
    }           
                                                                                                                                                                                                           
    cout << "      AGCM: init_tropopause_layers ended" << endl;                                                                                                                                          
}
/*
*
*/
// ============================================================================
// Atmosphere Model Initialization - Improved Version
// ============================================================================

#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>

// ============================================================================
// Physical Constants (should ideally be in a separate constants header)
// ============================================================================
namespace AtmosphereConstants {
    // Temperature constants
    constexpr double BETA_COSMO = 44.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 42.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 38.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 35.0;              // [K] COSMO parameter
//    constexpr double BETA_COSMO = 30.0;              // [K] COSMO parameter
    constexpr double MIN_SAFE_TEMP = 100.0;          // [K] Minimum safe temperature
    
    // Humidity thresholds
    constexpr double MIN_HUMIDITY = 0.0;             // [%]
    constexpr double MAX_HUMIDITY = 100.0;           // [%]
    
    // Water/air ratio factors
    constexpr double MIN_WATER_FACTOR = 0.5;         // Minimum total water factor
}

using namespace AtmosphereConstants;

// ============================================================================
// Helper Function: Project temperature to sea level (algebraic COSMO inversion)
// ============================================================================
/**
 * Inverts the COSMO temperature profile T(h) = sqrt(T0^2 - 2*beta*g*h/R)
 * exactly, so that feeding T0 back into the vertical reconstruction reproduces
 * t_u_init at h_mt without any round-trip error.
 *
 * T0 = sqrt(t_u_init^2 + 2*beta*g*h_mt / R_Air)
 *
 * @param t_u_init  Surface temperature at mountain height h_mt [K]
 * @param h_mt      Mountain height [m]
 * @return pair<T0_sea_level [K], p0_sea_level [hPa]>
 */
std::pair<double, double> project_to_sea_level(
    double t_u_init,
    double h_mt,
    double beta,
    double R_Air,
    double r_air,
    double g)
{
    double T0 = sqrt(t_u_init * t_u_init + (2.0 * beta * g * h_mt) / R_Air);
    double p0 = 1e-2 * (r_air * R_Air * T0);
    return {T0, p0};
}

// ============================================================================
// Main Initialization Function
// ============================================================================
void cAtmosphereModel::initTemperatureData(int Ma) {
    std::cout << "\n\n\n      AGCM: initTemperatureData" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // ATHAD: prescribed Hadean surface temperature
    // ========================================================================
    //
    // The inherited Earth code built t.x[0] from the NASA surface field (Ma == 0) or a
    // Scotese pole→equator parabola (Ma > 0), with EarthByte reconstruction corrections
    // between slices. None of that exists for 4.4 Ga, so ATHAD prescribes the surface
    // temperature directly from two parameters and keeps the parabola SHAPE only:
    //
    //     T_s(j) = (t_surf_pole - t_surf_equator) * parabola(ratio) + t_surf_pole
    //
    // with parabola(x) = x² - 2x, ratio = j / ((jm-1)/2), so ratio = 1 at the equator
    // (parabola = -1, giving t_surf_equator) and ratio = 0 or 2 at the poles (parabola = 0,
    // giving t_surf_pole). This is the identical algebra the paleo branch used, which keeps
    // the correspondence with ATOM_Precipitation readable.
    //
    // t.x[0] is stored non-dimensional (T / t_0), the convention Step 8 reads back.
    //
    // ASSUMPTION, not a result: a 250 bar steam atmosphere is optically thick enough that
    // the equator-pole surface contrast should be small, and 1500/1450 K is a guess at it.
    // See CLAUDE.md, "Stated assumptions".

    const double inv_t0    = 1.0 / t_0;
    const double d_j_half  = 0.5 * (jm - 1);

    double t_surf_sum = 0.0;

    #pragma omp parallel for collapse(2) reduction(+: t_surf_sum)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            const double ratio = (double)j / d_j_half;
            const double T_s   = (t_surf_pole - t_surf_equator) * AtomUtils::parabola(ratio)
                               + t_surf_pole;                            // [K]

            t.x[0][j][k] = T_s * inv_t0;                                 // non-dimensional
            t_surf_sum  += T_s;
        }
    }

    // Diagnostic mean, in °C to match the units the reporting code expects.
    t_global_mean = t_surf_sum / (double)(jm * km) - t_0;

    std::cout.precision(4);
    std::cout << "\n       ATHAD: prescribed Hadean surface temperature"
              << "\n       equator ......................................... t_surf_equator   = "
              << t_surf_equator << " K"
              << "\n       pole ............................................ t_surf_pole      = "
              << t_surf_pole << " K"
              << "\n       area mean ....................................... t_global_mean    = "
              << t_global_mean + t_0 << " K\n\n";

    // ========================================================================
    // Step 8: Vertical Temperature Profile & Potential Temperature
    // ========================================================================
    const double R_mix = m_comp.R_mix;
    const double R_bg  = m_comp.R_bg;
    const double M_bg  = m_comp.M_bg;

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            int i_mount = i_topography[j][k];
            double t_u_init = t.x[0][j][k] * t_0;                       // [K]

            // ATHAD: no land, so no land-sea thermal contrast and no elevation to project
            // through. i_mount is 0 in every column, so the surface IS the reference level
            // and t_u_init stands as prescribed by initTemperatureData.
            //
            // Surface pressure uses R_mix, the COLUMN reference gas constant, because r_air
            // was calibrated as p/(R_mix*T) — using R_Air (the non-condensable background,
            // 317.3) here instead gives 204 bar rather than the intended 250.
            p_stat.x[0][j][k] = 1e-2 * (r_air * R_mix * t_u_init);
            t.x[0][j][k]      = t_u_init;

            temp_pot.y[j][k]     = t.x[0][j][k];
            temp_reconst.y[j][k] = t.x[0][j][k];

            // ================================================================
            // Surface Humidity Calculation (Magnus Formula)
            // ================================================================
            const double t_u_0 = t_u_init;
            const double p_u_0 = p_stat.x[0][j][k];
            
            // Magnus coefficients depend on phase (water vs ice)
            const double E_Satur = SaturationH2O::saturationPressureAuto(t_u_0);
            const double e_curr  = c.x[0][j][k] * p_u_0 / (c.x[0][j][k] + ep);

            relative_humidity.y[j][k] = std::clamp(
                (e_curr / E_Satur) * 100.0,
                MIN_HUMIDITY,
                MAX_HUMIDITY
            );

            // ================================================================
            // Vertical Profile: Temperature, Pressure, Density
            // ================================================================
            // ATHAD: dry adiabat + isothermal skin, integrated hydrostatically.
            // Identical scheme to ThermoAtm::densities() — see the long note there for why
            // the inherited COSMO sqrt profile was replaced rather than re-tuned.
            double T_prev = std::max(MIN_SAFE_TEMP, t.x[0][j][k]);
            double p_prev = p_stat.x[0][j][k];

            for (int i = 0; i < im; i++) {
                const double q_v   = c.x[i][j][k];
                const double q_c   = co2.x[i][j][k];
                const double R_loc = AtmMixture::R_of(q_v, q_c, R_bg);

                double T_curr, p_val;
                if (i == 0) {
                    T_curr = T_prev;
                    p_val  = p_prev;
                } else {
                    const double dz     = get_layer_height(i) - get_layer_height(i-1);
                    const double cp_loc = AtmMixture::cp_of(q_v, q_c, T_prev, M_bg);
                    const double T_ad   = T_prev - (g / cp_loc) * dz;
                    T_curr = std::max(t_skin, T_ad);
                    const double T_mean = 0.5 * (T_prev + T_curr);
                    p_val  = p_prev * exp(-g * dz / (R_loc * T_mean));
                }

                t.x[i][j][k]      = T_curr;
                p_stat.x[i][j][k] = p_val;

                const double total_water_factor =
                    std::max(MIN_WATER_FACTOR, 1.0 - (cloud.x[i][j][k] + ice.x[i][j][k]));
                const double R_dloc = AtmMixture::R_of(0.0, q_c, R_bg);

                r_dry.x[i][j][k]   = p_val * 100.0 / (R_dloc * T_curr);
                r_humid.x[i][j][k] = p_val * 100.0 / (R_loc * T_curr * total_water_factor);

                T_prev = T_curr;
                p_prev = p_val;
            }

            temp_landscape.y[j][k]   = t.x[i_mount][j][k] - t_0;        // [°C]
            p_stat_landscape.y[j][k] = p_stat.x[i_mount][j][k];
        }
    }

    // ========================================================================
    // Step 9: Convert to Non-Dimensional and Apply Boundary Conditions
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            temp_reconst.y[j][k] = t.x[0][j][k] - t_0;  // [°C]

            // Convert all vertical levels to non-dimensional
            for (int i = 0; i < im; i++) {
                t.x[i][j][k] *= inv_t0;
            }

            // Apply topography boundary condition
            int i_mount = i_topography[j][k];
            if (i_mount >= 0 && i_mount < im) {
                t.x[0][j][k] = t.x[i_mount][j][k];
                p_stat.x[0][j][k] = p_stat.x[i_mount][j][k];
                r_dry.x[0][j][k] = r_dry.x[i_mount][j][k];
                r_humid.x[0][j][k] = r_humid.x[i_mount][j][k];
            }
        }
    }

    // ========================================================================
    // Step 9b (OPTION A): smooth the topography-draped surface temperature.
    // The DEM enters the lapse through the INTEGER level i_topography, so the draped
    // terrain-surface temperature has grid-scale steps (~one layer ~400 m ~2 K) at
    // cliffs/coasts -- exactly where this model's 2dx coastal/pressure modes historically
    // ignite. Apply a few light 1-2-1 passes (phi periodic, poles fixed) to the 2D surface
    // field and write it back into the PROGNOSTIC surface cell t.x[i_mount]; bcSolidGround
    // copies i_mount->0 every step, so smoothing i=0 alone would be undone on iteration 1.
    // Paleo only -- the modern NASA field and ocean cells (i_mount=0) are left untouched.
    if (*get_current_time() != 0) {
        const int n_passes = 4;
        std::vector<double> Tsurf(jm * km), Ttmp(jm * km);
        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++)
                Tsurf[j * km + k] = t.x[i_topography[j][k]][j][k];

        for (int p = 0; p < n_passes; p++) {
            for (int j = 0; j < jm; j++)                                // phi (k) pass, periodic
                for (int k = 0; k < km; k++) {
                    int km1 = (k - 1 + km) % km, kp1 = (k + 1) % km;
                    Ttmp[j * km + k] = 0.25 * Tsurf[j * km + km1]
                                     + 0.50 * Tsurf[j * km + k]
                                     + 0.25 * Tsurf[j * km + kp1];
                }
            for (int k = 0; k < km; k++) {                              // theta (j) pass, poles held
                Tsurf[k]              = Ttmp[k];
                Tsurf[(jm - 1) * km + k] = Ttmp[(jm - 1) * km + k];
                for (int j = 1; j < jm - 1; j++)
                    Tsurf[j * km + k] = 0.25 * Ttmp[(j - 1) * km + k]
                                      + 0.50 * Ttmp[j * km + k]
                                      + 0.25 * Ttmp[(j + 1) * km + k];
            }
        }

        for (int j = 0; j < jm; j++)
            for (int k = 0; k < km; k++) {
                int i_mount = i_topography[j][k];
                t.x[i_mount][j][k]     = Tsurf[j * km + k];             // prognostic surface cell
                t.x[0][j][k]           = Tsurf[j * km + k];             // sea-level reference layer
                temp_landscape.y[j][k] = Tsurf[j * km + k] * t_0 - t_0; // keep the diagnostic in sync [degC]
            }
    }

/*
    // ========================================================================
    // Post-processing: smooth t, p_stat, r_dry, r_humid horizontally                                        
    // ========================================================================
    {
        // Precompute Gaussian weights for the stencil (constant for all levels/cells)
        const int R = 2;
        const int D = 2 * R + 1;

        std::vector<double> gauss_w(D * D);

        for (int dj = -R; dj <= R; dj++)
            for (int dk = -R; dk <= R; dk++)

              gauss_w[(dj + R) * D + (dk + R)] =
                    std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));


        #pragma omp parallel
        {
            std::vector<double> tmp(jm * km);                           // thread-local scratch: one slice per thread

          // Lambda captures thread-local tmp — safe to call from any thread
            auto smoothLevel = [&](auto& field, int i, auto pred) {
                // Snapshot level i before writing (read-then-write on same array)
                for (int j = 0; j < jm; j++)
                    for (int k = 0; k < km; k++)
                        tmp[j * km + k] = field.x[i][j][k];

                for (int j = 0; j < jm; j++) {
                    for (int k = 0; k < km; k++) {
                        if (!pred(j, k)) continue;

                        double sum = 0.0;
                        double cnt = 0;

                        for (int dj = -R; dj <= R; dj++) {
                            const int jj = std::clamp(j + dj, 0, jm - 1);

                            for (int dk = -R; dk <= R; dk++) {
                                const int kk = (k + dk + km) % km;

                                if (i >= i_topography[jj][kk]) {
//                                    const double w = std::exp(-0.5 * (dj*dj + dk*dk) / double(R*R));   // if not use Gaussian weights
                                    const double w = gauss_w[(dj + R) * D + (dk + R)];
                                    sum += w * tmp[jj * km + kk];
                                    cnt += w;
                                }
                            }
                        }
                        field.x[i][j][k] = cnt > 0.0 ? sum / cnt : 0.0;
                    }
                }
            };

            // --- t and p_stat: above topography only ---
            // Each i-level is independent: no data dependency between levels
            #pragma omp for schedule(dynamic, 4)
            for (int i = 1; i < im - 1; i++) {
                auto above_topo = [&](int j, int k) {
                    return i >= i_topography[j][k];
                };

                smoothLevel(t,      i, above_topo);
                smoothLevel(p_stat, i, above_topo);
            }
            // implicit barrier: all t/p_stat levels done before r_dry/r_humid start

            // --- r_dry and r_humid: ocean cells only (land stays 0) ---
            #pragma omp for schedule(dynamic, 4)
            for (int i = 1; i < im - 1; i++) {
                auto ocean_above_topo = [&](int j, int k) {
                    return i >= i_topography[j][k] && !is_land(h, i, j, k);
                };

                smoothLevel(r_dry,   i, ocean_above_topo);
                smoothLevel(r_humid, i, ocean_above_topo);
            }
        } // end omp parallel — implicit barrier before BC loop

        // Re-apply i=0 boundary conditions from smoothed i_mount values
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const int im0 = i_topography[j][k];

                if (im0 >= 0 && im0 < im) {
                    t.x[0][j][k]       = t.x[im0][j][k];
                    p_stat.x[0][j][k]  = p_stat.x[im0][j][k];
                    r_dry.x[0][j][k]   = r_dry.x[im0][j][k];
                    r_humid.x[0][j][k] = r_humid.x[im0][j][k];
                }

                temp_reconst.y[j][k]     = t.x[0][j][k] * t_0 - t_0;
                temp_landscape.y[j][k]   = t.x[im0][j][k] * t_0 - t_0;
                p_stat_landscape.y[j][k] = p_stat.x[im0][j][k];
            }
        }
    }
*/


    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "      Initialization completed in " << duration.count() << " ms\n";
}
/*
*  ATHAD: the Scotese et al. (2021) paleo-temperature curve loaders and the
*  interpolating accessor get_temperatures_from_curve() lived here. They are gone:
*  the curves stop at ~540 Ma and nothing reaches 4.4 Ga, so ATHAD prescribes its
*  surface temperature from t_surf_equator / t_surf_pole instead (initTemperatureData).
*
*  The accessor was also unsafe as written — it dereferenced m.begin() and decremented
*  m.end() BEFORE its own m.size() < 2 guard, so calling it on an empty map was
*  undefined behaviour rather than the intended NAN return. With the curves removed
*  that path was live on every printDataAtm() call.
*/
/*
*
*/
// ATHAD: the surface is entirely water, so this only confirms the invariant.
//
// The Earth version also reported the per-point CO2 budget added by the ocean and land
// surfaces and removed by vegetation. The Hadean has no vegetation, no land and no
// carbonate ocean sink, so those terms are gone with their parameters. Note this
// function divided by h_land to form the ocean/land ratio — with no land that is a
// division by zero, which is the other reason it could not be left as it was.
void cAtmosphereModel::LandOceanFraction(){

    cout << endl << endl << endl << "      AGCM: LandOceanFraction" << endl;

    int h_land = 0;
    for(int j = 0; j < jm; j++)
        for(int k = 0; k < km; k++)
            if(is_land(h, 0, j, k)) h_land++;

    const int h_point_max = jm * km;

    cout.precision(3);
    cout << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.')
        << "      total number of surface points " << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << h_point_max << endl << setiosflags(ios::left) << setw(50)
        << setfill('.') << "      number of points on the water surface " << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << h_point_max - h_land << endl << setiosflags(ios::left) << setw(50)
        << setfill('.') << "      number of points on the land surface " << " = "
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ')
        << h_land << endl << endl;

    // Invariant 1 in CLAUDE.md: there is no topography, so is_land() must be false
    // everywhere. If that ever stops holding, every land branch in the RHS, the BCs,
    // the turbulence and the ice schemes silently comes back to life.
    if(h_land != 0){
        cout << "      ERROR: " << h_land << " land points found on a surface that must be "
             << "entirely water. ATHAD prescribes no topography." << endl;
        throw std::logic_error("   ATHAD: land points found on the flat Hadean surface");
    }
    cout << "      flat Hadean surface confirmed: 100% water, no land points" << endl << endl;

    cout << "      AGCM: LandOceanFraction ended" << endl;
}
/*
*
*/
void cAtmosphereModel::initWaterWapour() {
    std::cout << "\n\n\n      AGCM: initWaterWapour" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();


    // ========================================================================
    // Main Computation Loop: water vapour field
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            int i_mount = i_topography[j][k];

            // ATHAD: water vapour is WELL MIXED at its composition mass fraction.
            //
            // The Earth version set c to a fraction of the local saturation mixing ratio,
            // because on Earth water vapour is a condensable trace whose abundance IS set
            // by saturation. Here it is 67 % of the atmosphere's mass and, below the
            // condensation level, supercritical — there is no saturation to be a fraction
            // of. Evaluating the Magnus formula at 1500 K returns E_sat ~ 1.2e7 hPa, far
            // above the 250 bar total pressure, so the q_sat branch collapsed to its
            // fallback and the surface scheme then drove c to 20.8 — a mass fraction twenty
            // times larger than all the mass present.
            //
            // So c is initialised exactly as co2 is: uniform at the configured value. Where
            // the column does rise above the condensation level (T < 647 K, the top ~50 km),
            // the saturation adjustment will draw it down — once Phase 4 gives it a
            // saturation curve that is valid there.
            (void)i_mount;
            for (int i = 0; i < im; i++) {
                c.x[i][j][k]     = c_0;
                cloud.x[i][j][k] = 0.0;
            }
        }
    }

    // ========================================================================
    // Surface Boundary Condition (copy topography values to surface)
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
            if (i_mount >= 0 && i_mount < im && is_land(h, 0, j, k)) {
                c.x[0][j][k] = c.x[i_mount][j][k];
                cloud.x[0][j][k] = cloud.x[i_mount][j][k];
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initWaterWapour\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: initWaterWapour ended" << std::endl;
}
/*
*
*/
void cAtmosphereModel::initCloudIce() {
    std::cout << "\n\n\n      AGCM: initCloudIce" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double alfa_s    = 1.5;
    const double Hu_cr_max = 1.0;
    const double Hu_cr_mid = 0.8;
    const double Hu_diff   = Hu_cr_max - Hu_cr_mid;
//    const double det_T_0   = t_0 - 3.0;
    const double det_T_0   = t_0;
    // Parabola H_crit(p): roots at x=0 and x=1, minimum Hu_cr_mid at x = x_mid.
    // H_crit = Hu_cr_max - Hu_curv * x * (1 - x)
    //
    // x IS A FRACTION OF THE LOCAL SURFACE PRESSURE, not an absolute pressure. The inherited
    // form was x = p/1000 hPa with the minimum at 550 hPa — Earth's surface pressure and
    // Earth's mid-troposphere. Both are Earth constants, and at 250 bar they do something
    // much worse than going inert:
    //
    //   levels 0-35 (0-201 km)   p > 1000 hPa, so x > 1, the parabola goes NEGATIVE and
    //                            H_crit is saved only by the clamp below -> 1.0 exactly
    //   level 36    (218 km)     p = 439 hPa, x = 0.44 -> H_crit = 0.8010, the parabola's
    //                            MINIMUM, landing inside the condensable band (T < 647 K
    //                            only above level 34) and on the level where tau_above
    //                            crosses 1 -- the photosphere (README items 39, 41)
    //   levels 37-40             H_crit recovers 0.917 -> 0.9998
    //
    // So Earth's mid-troposphere cloud-threshold minimum was being applied at 218 km, making
    // cloud 20 % easier to form on precisely the level that sets the outgoing flux. That is
    // why initCloudIce was measured as the only live path from the grid to the OLR (item 41):
    // max cloud water moved 37.3 -> 49.9 g/kg with a change of grid and moist physics off.
    //
    // Rescaled, x = p/p_surface, the parabola's minimum moves to ~36 km — deep in the
    // supercritical column where nothing can condense — and H_crit >= 0.9926 everywhere
    // condensation IS possible. The Earth profile becomes inert here BY DERIVATION rather
    // than by fiat, which is the honest answer: a critical-relative-humidity fit to Earth's
    // cloud climatology has no analogue in a 250 bar steam shell whose condensable region is
    // a thin band at 185-300 km.
    //
    // EXACTLY EQUIVALENT ON EARTH, where p_surface = 1000 hPa, so this is portable upstream.
    // ATM_HCRIT_ABS=1 restores the inherited absolute-pressure form for A/B.
    const double x_mid   = 0.55;                                        // 550/1000 hPa on Earth
    const double Hu_curv = Hu_diff / (x_mid * (1.0 - x_mid));           // ~0.808

    static const bool hcrit_abs = [](){
        const char* e = getenv("ATM_HCRIT_ABS"); return e && atoi(e) != 0; }();
    if(hcrit_abs){
        std::cout << "      AGCM: initCloudIce - ATM_HCRIT_ABS set, H_crit on ABSOLUTE"
                  << " pressure (inherited Earth form, minimum at 550 hPa)" << std::endl;
    }

    // ========================================================================
    // Pass 1: cloud_max[i] — parallel over i, sum thread-local
    // ========================================================================
    cloud_max = std::vector<double>(im, 0.0);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < im; i++) {
        double sum = 0.0;
        for (int j = 0; j < jm; j++) {
            for (int k = 0; k < km; k++) {
                const double t_u = t.x[i][j][k] * t_0;
                const double p_u = p_stat.x[i][j][k];

                // ATHAD: nothing condenses above the critical point, and nothing condenses
                // where the vapour is superheated either.
                //
                // The critical-point guard below was the first half of this fix. The second
                // half is the conversion itself: q_sat = ep*E_sat/(p_u - E_sat) goes NEGATIVE
                // whenever E_sat exceeds the local pressure, and (c - H_crit*q_sat) then
                // manufactures cloud out of the subtraction of a negative number. The guard
                // stopped that above 647 K but not in the SUBcritical band between ~373 K and
                // the critical point, where p_sat still exceeds the local pressure through the
                // whole upper column — 140 to 200 km here. The exact form returns 1 there,
                // which is what "no condensation is possible" actually means.
                if (t_u >= AtmMixture::T_CRIT_H2O) continue;

                const double E_sat = (t_u >= t_0)
                    ? SaturationH2O::saturationPressure(t_u)
                    : SaturationH2O::sublimationPressure(t_u);
                const double q_sat = SaturationH2O::saturationMassFraction(
                    E_sat, p_u, AtmMixture::M_nonwater(c.x[i][j][k], co2.x[i][j][k],
                                                       m_comp.M_bg));

//                sum += std::max(0.0, c.x[i][j][k] - q_sat); }
//                sum += std::max(0.0, c.x[i][j][k] - 0.84 * q_sat); }
                sum += std::max(0.0, c.x[i][j][k] - 0.74 * q_sat); }
        }
        cloud_max[i] = sum / ((jm-1) * (km-1));
        if (cloud_max[i] <= 1e-4)  cloud_max[i] = 1e-4;
    }

    // ========================================================================
    // Precompute per-level quantities (sequential: im is small)
    // ========================================================================
    std::vector<double> step(im, 0.0);
    std::vector<double> dt_dim(im, 0.0);
    std::vector<double> alfa_over_cmax(im, 0.0);
    std::vector<double> two_step(im, 0.0);

    for (int i = 0; i < im - 1; i++) {
        step[i]           = get_layer_height(i+1) - get_layer_height(i);
        dt_dim[i]         = step[i] / 1.6;
        alfa_over_cmax[i] = alfa_s / cloud_max[i];
        two_step[i]       = 2.0 * step[i];
    }
    step[im-1]           = step[im-2];
    dt_dim[im-1]         = dt_dim[im-2];
    alfa_over_cmax[im-1] = alfa_over_cmax[im-2];
    two_step[im-1]       = two_step[im-2];

    // ========================================================================
    // Pass 2: cloud, ice, graupel fields — collapse(j, k), i inner
    // ========================================================================
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            // Normaliser for the H_crit parabola, once per column. p_stat.x[0] is the surface
            // and is p_0 everywhere (no topography, and item 19 pinned it), but read it rather
            // than assume it so a future column-varying surface pressure stays correct.
            const double p_surf_jk  = p_stat.x[0][j][k];
            const double inv_p_norm = hcrit_abs
                                    ? (1.0 / 1000.0)                    // inherited Earth form
                                    : ((p_surf_jk > 0.0) ? 1.0 / p_surf_jk : 1.0 / 1000.0);

            for (int i = 0; i < im; i++) {
                const double t_u = t.x[i][j][k] * t_0;
                const double p_u = p_stat.x[i][j][k];

                // ATHAD: nothing condenses above the critical point, and nothing condenses
                // where the vapour is superheated either.
                //
                // The critical-point guard below was the first half of this fix. The second
                // half is the conversion itself: q_sat = ep*E_sat/(p_u - E_sat) goes NEGATIVE
                // whenever E_sat exceeds the local pressure, and (c - H_crit*q_sat) then
                // manufactures cloud out of the subtraction of a negative number. The guard
                // stopped that above 647 K but not in the SUBcritical band between ~373 K and
                // the critical point, where p_sat still exceeds the local pressure through the
                // whole upper column — 140 to 200 km here. The exact form returns 1 there,
                // which is what "no condensation is possible" actually means.
                if (t_u >= AtmMixture::T_CRIT_H2O) {
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k]   = 0.0;
                    gr.x[i][j][k]    = 0.0;
                    continue;
                }

                const double E_sat = (t_u >= t_0)
                    ? SaturationH2O::saturationPressure(t_u)
                    : SaturationH2O::sublimationPressure(t_u);
                const double q_sat = SaturationH2O::saturationMassFraction(
                    E_sat, p_u, AtmMixture::M_nonwater(c.x[i][j][k], co2.x[i][j][k],
                                                       m_comp.M_bg));

                const double x_norm = p_u * inv_p_norm;
                double H_crit = Hu_cr_max - Hu_curv * x_norm * (1.0 - x_norm);
                if (H_crit > 1.0)  H_crit = 1.0;

                const double del_q_ls = std::max(0.0, c.x[i][j][k] - H_crit * q_sat);
                const double cloud_ls = cloud_max[i] * (1.0 - exp(-alfa_over_cmax[i] * del_q_ls));

                double cloud_conv = 0.0;
                if (P_rain.x[i][j][k] > 0.0 && i < im - 1) {
                    const double del_q_conv = std::max(0.0,
                        (P_rain.x[i+1][j][k] - P_rain.x[i][j][k])
                        / (two_step[i] * r_humid.x[i][j][k]) * dt_dim[i]);
                    cloud_conv = cloud_max[i] * (1.0 - exp(-alfa_over_cmax[i] * del_q_conv));
                }

                double cloud_val = cloud_ls + cloud_conv;
                if (is_land(h, i, j, k))  cloud_val = 0.0;
                cloud.x[i][j][k] = cloud_val;
 
                double h_T = 0.0;
                if (t_u < t_0) {
                    const double ratio = (t_u - t_0) / det_T_0;
                    h_T = 1.0 - exp(-0.5 * ratio * ratio);
                }

                ice.x[i][j][k] = cloud_val * h_T;
                gr.x[i][j][k]  = 0.1 * cloud_val * h_T;
 
                if (t_u <= t_00) {
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k]   = 0.0;
                    gr.x[i][j][k]    = 0.0;
                }
            }  // i
        }  // k
    }  // j

    // ========================================================================
    // Surface Boundary Condition
    // ========================================================================
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < jm; j++) {
        for (int k = 0; k < km; k++) {
            const int i_mount = i_topography[j][k];
//            if (i_mount >= 0 && i_mount < im && is_land(h, 0, j, k)) {
            if (i_mount >= 0 && i_mount < im && is_land(h, i_mount, j, k)) {
                cloud.x[0][j][k] = cloud.x[i_mount][j][k];
                ice.x[0][j][k]   = ice.x[i_mount][j][k];
                gr.x[0][j][k]    = gr.x[i_mount][j][k];
            }
            for (int i = i_mount-1; i >= 0; i--) {
                if (is_land(h, i, j, k)) {
                    cloud.x[i][j][k] = 0.0;
                    ice.x[i][j][k]   = 0.0;
                    gr.x[i][j][k]    = 0.0;
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" Time measured: %.3f seconds for initCloudIce\n", elapsed.count() * 1e-9);

    std::cout << "      AGCM: initCloudIce ended" << std::endl;
}
/*
*
*/

/*
* initBalancedState — ported from ASTIM cf43bfd ("Stop the prescribed jet from spinning
* down"), which found two independent causes of a prescribed circulation decaying. The
* second, the radial Shapiro pass draining vertical shear, was A/B-tested here and
* ACQUITTED: turning it fully off (ATM_RADIAL_SHAPIRO_STRENGTH=0) changes Psi_max by
* 0.16 % at 20 iterations and the mid-latitude cell still dies at the same iteration.
* That annotation is true where it was written -- Earth, sharp jets on a 16 km shell --
* and false on 41 levels spread over 300 km. This is the FIRST cause, which ATHAD has
* no answer for at all.
*
* A prescribed circulation is not a solution of anything on its own. VelocityInitializer
* imposes an analytic Hadley/Ferrel profile onto a field that carries no pressure
* structure to support it, so from iteration 0 the Coriolis torque on the imposed wind is
* unopposed and the flow accelerates instead of the cell equilibrating. Measured here:
* at 45 deg the surface streamfunction crosses zero at iteration ~5 and then grows
* LINEARLY at ~1335 per iteration -- constant dv/dt, the free-acceleration signature --
* burying the two-branch cell structure within five iterations.
*
* THE BALANCE, TAKEN FROM THIS MODEL'S OWN theta-MOMENTUM EQUATION rather than an
* idealised one, which is the part of ASTIM's approach worth keeping. With u = v = 0 and
* d/dphi = 0, RHS_Atm_Turb.cpp writes
*
*     rhs_v = -dpdthe*inv_rm - transport_v + coriolis*force_nd*coriolis_the + ...
*     transport_v  += (u*v - w^2*cotanthe) * inv_rm          [line 826]
*     coriolis_the  = 2*costhe*w                             [line 394]
*
* and requiring rhs_v = 0 gives
*
*     dp_dyn/dthe = w^2*cotanthe + 2*force_nd*costhe*w*rm
*
* with force_nd = omega*metricShellLength()/u_0 and the sin(theta) >= 0.55 metric floor,
* i.e. the same coefficients and the same polar clamp the RHS uses. ASTIM needed a
* predictor in its trapezoid because its integrand depended on the density it was
* solving for; ATHAD's depends only on w and the geometry, so a plain trapezoid is exact.
*
* WHY IT GOES INTO p_dyn, AND WHY IT SURVIVES. p_stat appears nowhere in the momentum
* equations (only in the Held-Suarez sigma), so the ENTIRE meridional pressure-gradient
* force in this model comes from p_dyn -- a 50 K equator-to-pole contrast over a 250 bar
* column exerts none of it directly. p_dyn is the projection variable, and
* project_initial_velocity() zeroes it, so this must run AFTER that call; the time loop's
* run() then relaxes p_dyn in place rather than replacing it. This is also why item 18
* found pgf stuck at 1 % of Coriolis after 400 iterations: a geostrophically balanced
* pressure is very nearly divergence-free, so it sits in the null space of what the
* projection solves each step, and the flow has to build it the slow way through the
* buoyancy-driven divergence. Supplying it at t = 0 is the point.
*
* The integration constant is fixed by removing the sin(theta)-weighted mean at each
* level, so the balance redistributes pressure in latitude without adding any.
*
* ATM_BALANCED_INIT scales the whole perturbation. DEFAULT 0.0 = OFF and bit-identical,
* per this repo's convention for a change that has not yet been measured over a long run;
* 1.0 is full balance, and intermediate values exist because ASTIM found the response
* strongly nonlinear (98.3 % retention at 0.0 against 92.7 % at 1.0 for its filter knob).
*
* ============================================================================
* THE RADIAL HALF, AND WHY THE THETA-ONLY VERSION ABOVE IS NOT A BALANCE (item 28).
*
* Everything above balances ONE component. p_dyn is one scalar field and it appears in
* all three momentum equations, so fixing rhs_v = 0 fixes dp/dtheta and says nothing
* about dp/dr -- and the radial equation then gets whatever is left. It is not small.
* Measured over 200 iterations at kappa_H2O = 0.010, balanced against the unbalanced
* control, both with the moist physics from iteration 0:
*
*     max |u|        0.114 m/s (iter 0)  ->  11.17 m/s (iter 200)     balanced
*                    0.114 m/s           ->   0.095 m/s               unbalanced
*     u at the equator, 95 km   +0.061 -> -6.35 m/s: the tropics REVERSE to sinking
*
* The aspect ratio says the vertical velocity should be ~(300 km/10 000 km)*v ~ 0.09 m/s,
* which is what the control has. The balanced run is ~100x over it, with descent over the
* 1500 K equator and ascent at the poles. Psi_max -- the number item 27 flipped the
* default on -- is built from v alone and reports 157 555 against 153 301 throughout:
* item 26's lesson recurring, an instrument that cannot see the component going wrong.
*
* THE MECHANISM. p_bal is integrated level by level and its amplitude follows w^2, which
* grows strongly upward, so the field varies with r even though each level's mean is
* removed. Measured at the equator: 6.9e6 hPa at 9 km -> 1.79e7 hPa at 55 km, i.e.
* -(1/rho)dp/dz of -630 to -1290 m/s^2. Nothing opposes it: with the DEFAULT switches
* the radial equation has no other term at u = v = 0. coriolis_rad carries the factor
* nontrad, and AtomUtils::coriolis_nontraditional() is false by default; the -(v^2+w^2)/r
* curvature term is inside if(AtomUtils::metric_curvature()), false by default; and the
* buoyancy term is multiplied by buoyancy_ramp, which is 0 at iteration 0. So
*
*     rhs_u = -dp_dyn/dr * exp_rm      and nothing else.
*
* WHICH ALSO CONVICTS THE THETA HALF. That same switch governs the w^2*cotanthe term the
* balance above is built on: with metric_curvature() off it is not in rhs_v either, so
* the theta balance is being computed against an equation the model does not solve. The
* balance was derived from the RHS as written, not as configured -- the family's constant
* trap in a new form, an Earth-tested code path read without checking whether it runs.
*
* THE FIX. There is no p(r,theta) satisfying both components exactly: that needs
* d(F_theta)/dr = d(F_r)/dtheta, and the curl of the force field is generally nonzero --
* the rotational part is what buoyancy balances in a real atmosphere, and here buoyancy is
* ramped off at iteration 0. So take the best available compromise instead of pretending
* one exists, minimising the acceleration the model is actually left holding:
*
*     J = SUM_ij [ (dp/dr*exp_rm - F_r)^2 + (dp/dtheta*inv_rm - F_theta)^2 ]
*
* in the RHS's own units -- both residuals ARE accelerations there, so this is the total
* spurious force. Its Euler-Lagrange equation is elliptic,
*
*     d/dr[ exp_rm^2 dp/dr ] + d/dtheta[ inv_rm^2 dp/dtheta ]
*         = d/dr[ exp_rm F_r ] + d/dtheta[ inv_rm F_theta ]
*
* solved here by SOR in conservative form, with the Neumann condition (dp/dr*exp_rm = F_r)
* at both walls falling out of dropping the boundary faces from both sides. F_r and
* F_theta are read from the RHS TERM BY TERM AND SWITCH BY SWITCH, so with the defaults
* F_r = 0 and the solve says what the radial equation says: p_dyn must not vary with r.
* Turn ATOM_METRIC_CURVATURE on and both halves gain their w^2 terms and the solve balances
* those instead -- the point of the formulation is that it does not have to be rewritten
* when the model is reconfigured.
*
* The gauge is now a single global constant, NOT a per-level mean: removing a different
* constant at each level is itself a radial gradient, which is the thing being fixed.
*
* ATM_BALANCED_MODE=1 restores the theta-only field for A/B (exactly as committed in item
* 27, curvature term unconditional, per-level gauge); 2 is this solve and is the default.
* Both print the residual accelerations they leave behind, so neither can be judged again
* by a diagnostic that only looks at one component.
*/
/*
* balancedStateSolve — the elliptic half of initBalancedState (README item 28).
*
* Given the non-pressure forces (F_r, F_theta) that the RHS applies at u = v = 0, find the
* p_dyn perturbation that leaves the least unbalanced acceleration in BOTH components:
*
*     minimise  J = SUM_ij [ (dp/dr*exp_rm - F_r)^2 + (dp/dtheta*inv_rm - F_theta)^2 ]
*     Euler-Lagrange:
*        d/dr[ exp_rm^2 dp/dr ] + d/dtheta[ inv_rm^2 dp/dtheta ]
*            = d/dr[ exp_rm F_r ] + d/dtheta[ inv_rm F_theta ]
*
* Conservative 5-point discretisation on the computational (i,j) grid, coefficients at the
* half points, boundary faces dropped from BOTH sides — which is exactly the Neumann
* condition dp/dr*exp_rm = F_r at the walls, i.e. the balance is asked for at the surface
* and the lid too rather than being clamped there.
*
* Solved by serial SOR. Serial ON PURPOSE: it is 41x181 = 7421 unknowns and costs
* milliseconds, and a deterministic sweep order keeps the initial state bit-identical run to
* run, which a red-black OpenMP sweep would not (README item 18 is about exactly this).
*
* Pure Neumann leaves the constant undetermined, so the mean is removed each sweep to stop
* it drifting, and once more at the end: that is the gauge, and it is GLOBAL. The theta-only
* path removed a separate mean at every level, and a per-level constant is a radial gradient.
*/
void cAtmosphereModel::balancedStateSolve(const std::vector<double>& F_r,
                                          const std::vector<double>& F_the,
                                          std::vector<double>& p) const
{
    const int    nij      = im * jm;
    const double inv_dr2  = 1.0 / (dr * dr);
    const double inv_dthe2 = 1.0 / (dthe * dthe);

    std::vector<double> e_half(im, 0.0);                      // exp_rm at i+1/2
    for(int i = 0; i < im - 1; i++)
        e_half[i] = 1.0 / (0.5 * (rad.z[i] + rad.z[i+1]) + 1.0);
    std::vector<double> inv_rm(im, 0.0);
    for(int i = 0; i < im; i++) inv_rm[i] = 1.0 / metricRadius(rad.z[i]);

    std::vector<double> cE(nij, 0.0), cW(nij, 0.0), cN(nij, 0.0), cS(nij, 0.0),
                        diag(nij, 0.0), src(nij, 0.0);

    for(int i = 0; i < im; i++){
        const double m2 = inv_rm[i] * inv_rm[i] * inv_dthe2;
        for(int j = 0; j < jm; j++){
            const int    id = i*jm + j;
            double s = 0.0;
            if(i < im - 1){
                cE[id] = e_half[i] * e_half[i] * inv_dr2;
                s += e_half[i] * 0.5 * (F_r[id] + F_r[id + jm]) / dr;
            }
            if(i > 0){
                cW[id] = e_half[i-1] * e_half[i-1] * inv_dr2;
                s -= e_half[i-1] * 0.5 * (F_r[id] + F_r[id - jm]) / dr;
            }
            if(j < jm - 1){
                cN[id] = m2;
                s += inv_rm[i] * 0.5 * (F_the[id] + F_the[id + 1]) / dthe;
            }
            if(j > 0){
                cS[id] = m2;
                s -= inv_rm[i] * 0.5 * (F_the[id] + F_the[id - 1]) / dthe;
            }
            diag[id] = cE[id] + cW[id] + cN[id] + cS[id];
            src [id] = s;
        }
    }

    // ALTERNATING-DIRECTION LINE RELAXATION, not point SOR. The operator is strongly anisotropic —
    // exp_rm^2/dr^2 runs ~55x the inv_rm^2/dthe^2 term, because a radial pressure gradient
    // buys ~10x the acceleration a meridional one does and the radial grid is finer. Point
    // relaxation crawls on that (20 000 sweeps and still moving); solving each radial column
    // exactly with Thomas and sweeping in theta converges in tens of sweeps, because the
    // stiff direction is no longer being iterated at all.
    const int    max_sweep = 5000;
    const double tol       = 1.0e-10;
    static const bool trace = (getenv("ATM_BALANCE_TRACE") != nullptr);
    double scale = 0.0;
    for(int id = 0; id < nij; id++) scale = std::max(scale, std::fabs(src[id]));
    if(!(scale > 0.0)){ std::fill(p.begin(), p.end(), 0.0); return; }

    const int nmax = std::max(im, jm);
    std::vector<double> a(nmax), b(nmax), c(nmax), d(nmax), cp(nmax), dp(nmax);

    // Thomas, in place on a/b/c/d, writing the solution back through `put`.
    // Over-relaxed: the line solves handle each direction exactly, but the error mode that
    // is smooth in BOTH is nearly in the null space of either pass and decays at ~0.999 per
    // sweep without this (measured: 5.5 -> 5.3e-3 over 5000 sweeps, still not converged).
    // 1.98 is the measured optimum on this grid: sweeps to a 1e-10 relative residual run
    // 5000+ (w=1.0, never gets there), 4265 (1.90), 1247 (1.98), 2441 (1.99), 5000+ (1.995).
    // The converged field is identical to 4 digits across all of them, as it must be — this
    // buys wall clock, not a different answer.
    static const double omega_line = [](){
        const char* e = getenv("ATM_BALANCE_OMEGA"); return e ? atof(e) : 1.98; }();
    auto thomas = [&](int n, auto&& get_old, auto&& put){
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        for(int m = 1; m < n; m++){
            const double q = b[m] - a[m] * cp[m-1];
            cp[m] = (std::fabs(q) > 1e-300) ? c[m] / q : 0.0;
            dp[m] = (std::fabs(q) > 1e-300) ? (d[m] - a[m] * dp[m-1]) / q : 0.0;
        }
        double dloc = 0.0, next = 0.0;
        for(int m = n - 1; m >= 0; m--){
            const double ps  = (m == n - 1) ? dp[m] : dp[m] - cp[m] * next;
            const double p0  = get_old(m);
            const double pn  = p0 + omega_line * (ps - p0);
            const double del = pn - p0;
            next = ps;                    // the sweep marches on the exact line solution
            put(m, pn);
            if(std::fabs(del) > dloc) dloc = std::fabs(del);
        }
        return dloc;
    };

    int    sweep = 0;
    double dmax  = 0.0;
    for(; sweep < max_sweep; sweep++){
        dmax = 0.0;

        // --- pass 2 of the alternating direction: theta lines, one per level. The radial
        // pass below kills the stiff direction; without this one the smooth theta mode
        // decays at ~1-(pi/jm)^2 per sweep and 5000 sweeps are not enough.
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                const int id = i*jm + j;
                a[j] = (j > 0)      ? -cS[id] : 0.0;
                c[j] = (j < jm - 1) ? -cN[id] : 0.0;
                b[j] = diag[id];
                d[j] = -src[id]
                     + (i < im-1 ? cE[id] * p[(i+1)*jm + j] : 0.0)
                     + (i > 0    ? cW[id] * p[(i-1)*jm + j] : 0.0);
                if(!(b[j] > 0.0)){ b[j] = 1.0; a[j] = c[j] = d[j] = 0.0; }
            }
            const double dl = thomas(jm,
                [&](int j){ return p[i*jm + j]; },
                [&](int j, double v){ p[i*jm + j] = v; });
            if(dl > dmax) dmax = dl;
        }

        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                const int id = i*jm + j;
                a[i] = (i > 0)      ? -cW[id] : 0.0;          // sub-diagonal
                c[i] = (i < im - 1) ? -cE[id] : 0.0;          // super-diagonal
                b[i] = diag[id];
                d[i] = -src[id]
                     + (j < jm-1 ? cN[id] * p[i*jm + j + 1] : 0.0)
                     + (j > 0    ? cS[id] * p[i*jm + j - 1] : 0.0);
                if(!(b[i] > 0.0)){ b[i] = 1.0; a[i] = c[i] = d[i] = 0.0; }
            }
            // The pure-Neumann column would be singular on its own (b = -(a+c) exactly when
            // there is no theta coupling); the diagonal carries the cN/cS terms, which make
            // it strictly dominant. Only an all-zero column can fail, and that is the guard.
            const double dl = thomas(im,
                [&](int i){ return p[i*jm + j]; },
                [&](int i, double v){ p[i*jm + j] = v; });
            if(dl > dmax) dmax = dl;
        }
        double mean = 0.0;
        for(int id = 0; id < nij; id++) mean += p[id];
        mean /= (double)nij;
        for(int id = 0; id < nij; id++) p[id] -= mean;

        // Stop on the RESIDUAL, not on the update. A pure-Neumann system relaxed with a
        // per-sweep mean removal keeps producing a constant-size update long after the
        // answer stops moving — measured here: the field was identical to four digits at
        // 100 sweeps and at 5000, while "last update" sat at 7e-4 and declared failure.
        // ||A p - src|| is the quantity that actually says whether it is solved.
        dmax = 0.0;
        for(int i = 0; i < im; i++){
            for(int j = 0; j < jm; j++){
                const int id = i*jm + j;
                if(!(diag[id] > 0.0)) continue;
                const double r = (i < im-1 ? cE[id] * p[id + jm] : 0.0)
                               + (i > 0    ? cW[id] * p[id - jm] : 0.0)
                               + (j < jm-1 ? cN[id] * p[id + 1]  : 0.0)
                               + (j > 0    ? cS[id] * p[id - 1]  : 0.0)
                               - diag[id] * p[id] - src[id];
                if(std::fabs(r) > dmax) dmax = std::fabs(r);
            }
        }
        if(trace && (sweep < 3 || sweep == 9 || sweep == 99 || sweep == 999
                     || sweep == 4999))
            std::cout << "        [balance trace] sweep " << (sweep + 1) << "  residual "
                      << std::scientific << std::setprecision(3) << dmax << std::fixed
                      << std::endl;
        if(dmax <= tol * scale) { sweep++; break; }
    }

    std::cout << "      balance solve: " << sweep << " line-relaxation sweeps, residual "
              << std::scientific << std::setprecision(3) << dmax << " against a source of "
              << scale << std::fixed
              << (sweep >= max_sweep ? "   <-- HIT THE SWEEP CAP, not converged" : "")
              << std::endl;
}
/*
*
*/
void cAtmosphereModel::initBalancedState(){
    // DEFAULT 1.0 = ON since README item 27. It was committed at 0.0 (bit-identical) and
    // flipped on a 200-iteration moist run at kappa_H2O = 0.010: Psi_max dips 0.34 % by
    // iteration 100 and returns to its starting value by 200, an oscillation instead of a
    // decay, while the unbalanced control lost 27 % of its tropical cell and grew a
    // mid-latitude drift to twice the cell's strength. ATM_BALANCED_INIT=0 forces the old
    // path for A/B. Intermediate values scale the perturbation.
    static const double strength = [](){
        const char* e = getenv("ATM_BALANCED_INIT"); return e ? atof(e) : 1.0; }();
    if(strength == 0.0) return;

    std::cout << std::endl << "      AGCM: initBalancedState begin ......................." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    // ATM_BALANCED_MODE — 2 (default) balances both momentum components at once; 1 is the
    // theta-only field of item 27, kept verbatim so the two can be run against each other.
    static const int mode = [](){
        const char* e = getenv("ATM_BALANCED_MODE"); return e ? atoi(e) : 2; }();

    const int j_eq = (jm - 1) / 2;
    const double force_nd = omega * metricShellLength() / u_0;

    // Count what the solver's clamp would truncate, read from the same accessor the solver
    // enforces so the two cannot drift apart. A balance the model then clips is not a
    // balance, and unlike a uniform scaling, clipping distorts the SHAPE. With the inherited
    // Earth ceiling of 10 this reported 94.7 % of columns over; item 27 raised it to 2000.
    double max_add = 0.0;
    long n_over = 0, n_tot = 0;
    const double ceiling_probe = pDynCeiling();

    // ---- the zonal-mean zonal wind the balance is built on, over fluid cells ----
    const int nij = im * jm;
    std::vector<double> wbar(nij, 0.0), F_r(nij, 0.0), F_the(nij, 0.0),
                        pbal2(nij, 0.0), padd(nij, 0.0);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            double ws = 0.0; int n = 0;
            for(int k = 0; k < km; k++){
                if(i < i_topography[j][k]) continue;          // fluid cells only
                ws += w.x[i][j][k]; n++;
            }
            wbar[i*jm + j] = (n > 0) ? ws / (double)n : 0.0;
        }
    }

    // ---- the forces the RHS actually applies at u = v = 0, d/dphi = 0 ----
    // Read from RHS_Atm_Turb.cpp term by term AND SWITCH BY SWITCH: coriolis_rad carries
    // nontrad, the curvature pair sits inside if(metric_curvature()), and both default off.
    // Signs follow rhs_x = ... - transport_x, so -transport_u = +(v^2+w^2)*inv_rm at v = 0.
    const bool   curvature = AtomUtils::metric_curvature();
    const double nontrad   = AtomUtils::coriolis_nontraditional() ? 1.0 : 0.0;

    for(int i = 0; i < im; i++){
        const double rmet   = metricRadius(rad.z[i]);
        const double inv_rm = 1.0 / rmet;
        for(int j = 0; j < jm; j++){
            double sinthe = sin(the.z[j]);
            if(sinthe < 0.55) sinthe = 0.55;                  // the RHS's metric floor
            const double costhe   = cos(the.z[j]);
            const double cotanthe = costhe / sinthe;
            const double wb       = wbar[i*jm + j];

            // No Coriolis flag here on purpose: RHS_Atm_Turb.cpp:375 hard-codes its own
            // local `coriolis = 1.0` and never reads the config member of that name, so
            // multiplying by the config value would balance a force the model does not apply.
            double fr  = force_nd * nontrad * 2.0 * sinthe * wb;
            double fth = force_nd * 2.0 * costhe * wb;
            if(curvature){
                fr  += wb * wb * inv_rm;
                fth += wb * wb * cotanthe * inv_rm;
            }
            F_r  [i*jm + j] = fr;
            F_the[i*jm + j] = fth;
        }
    }

    if(mode == 2) balancedStateSolve(F_r, F_the, pbal2);

    #pragma omp parallel for schedule(static) reduction(max: max_add) reduction(+: n_over, n_tot)
    for(int i = 0; i < im; i++){
        const double rm = rad.z[i];

        std::vector<double> pbal(jm, 0.0), f(jm, 0.0);

        // THE RADIUS HERE IS metricRadius(rm), NOT rad.z[i], and getting that wrong cost a
        // factor of ~20. The RHS multiplies BOTH the pressure gradient and the w^2*cotanthe
        // curvature term by inv_rm = 1/metricRadius(rm) — "the 1/r factors carry the
        // PLANETARY radius; exp_rm carries the stretched grid coordinate", as
        // RungeKutta_Atm_Turb.cpp puts it — while the Coriolis term carries no radius at all.
        // So for -dpdthe*inv_rm to cancel a Coriolis term that was never divided, dp/dthe has
        // to carry metricRadius(rm) = m_metric_r0 + (rm - rad.z[0]) ~ 21 in shell units, not
        // rm ~ 1. Measured with rm: pgf came out at 5.1 % of cor instead of 100 %. The
        // curvature term takes no factor, because it is divided by the same inv_rm the
        // pressure gradient is. This is the metric-disagreement defect that file warns about,
        // reproduced in the act of porting a fix.
        const double rmet = metricRadius(rm);
        double mean = 0.0;

        if(mode == 2){
            // The two-component solve has already produced the whole field, and its gauge
            // is one global constant (removed below), not a per-level mean.
            for(int j = 0; j < jm; j++) pbal[j] = pbal2[i*jm + j];
        }
        else{
            // Item 27's theta-only path, kept verbatim — including the unconditional
            // w^2*cotanthe term, which is what was measured, switch or no switch.
            for(int j = 0; j < jm; j++){
                const double the_j  = the.z[j];
                double sinthe = sin(the_j);
                if(sinthe < 0.55) sinthe = 0.55;              // the RHS's metric floor
                const double costhe   = cos(the_j);
                const double cotanthe = costhe / sinthe;
                const double wb       = wbar[i*jm + j];
                f[j] = wb * wb * cotanthe
                     + 2.0 * force_nd * costhe * wb * rmet;
            }

            for(int j = j_eq; j < jm - 1; j++)
                pbal[j+1] = pbal[j] + 0.5 * dthe * (f[j] + f[j+1]);
            for(int j = j_eq; j > 0; j--)
                pbal[j-1] = pbal[j] - 0.5 * dthe * (f[j] + f[j-1]);

            double wsum = 0.0, num = 0.0;
            for(int j = 0; j < jm; j++){
                const double aw = sin(the.z[j]);
                wsum += aw;
                num  += aw * pbal[j];
            }
            mean = (wsum > 0.0) ? num / wsum : 0.0;
        }

        for(int j = 0; j < jm; j++){
            const double add = strength * (pbal[j] - mean);
            padd[i*jm + j] = add;                             // kept for the residual report
            if(std::fabs(add) > max_add) max_add = std::fabs(add);
            n_tot++;
            if(std::fabs(add) > ceiling_probe) n_over++;
            for(int k = 0; k < km; k++)
                p_dyn.x[i][j][k] += add;
        }
    }

    std::cout << "      balancing pressure perturbation: mode " << mode
              << (mode == 2 ? " (two-component least squares)" : " (theta only, item 27)")
              << ", max |dp_dyn| = " << max_add
              << " (non-dim), strength = " << strength << std::endl;
    std::cout << "      above the p_dyn ceiling (" << ceiling_probe << "): " << n_over
              << " of " << n_tot << " (i,j) columns ("
              << (n_tot > 0 ? 100.0 * (double)n_over / (double)n_tot : 0.0)
              << " %) -- these would be clipped on the solver's first call" << std::endl;

    // WHAT THE BALANCE LEAVES BEHIND, IN BOTH COMPONENTS. rhs_u and rhs_v apply
    // -dp/dr*exp_rm and -dp/dtheta*inv_rm, so these residuals are the accelerations the
    // model is still holding at iteration 0, in its own units. Printing only the theta one
    // is how a field that drove the vertical wind to 11 m/s passed for a balance (item 28).
    {
        const double inv_2dr_l   = 1.0 / (2.0 * dr);
        const double inv_2dthe_l = 1.0 / (2.0 * dthe);
        double r2 = 0.0, t2 = 0.0, r2_0 = 0.0, t2_0 = 0.0, rmx = 0.0, tmx = 0.0;
        long n = 0;
        for(int i = 1; i < im - 1; i++){
            const double exp_rm = 1.0 / (rad.z[i] + 1.0);
            const double inv_rm = 1.0 / metricRadius(rad.z[i]);
            for(int j = 1; j < jm - 1; j++){
                const double dpdr   = (padd[(i+1)*jm + j] - padd[(i-1)*jm + j]) * inv_2dr_l;
                const double dpdthe = (padd[i*jm + j+1]   - padd[i*jm + j-1])   * inv_2dthe_l;
                const double res_r  = -dpdr   * exp_rm + F_r  [i*jm + j];
                const double res_t  = -dpdthe * inv_rm + F_the[i*jm + j];
                r2   += res_r * res_r;      t2   += res_t * res_t;
                r2_0 += F_r[i*jm+j] * F_r[i*jm+j];
                t2_0 += F_the[i*jm+j] * F_the[i*jm+j];
                if(std::fabs(res_r) > rmx) rmx = std::fabs(res_r);
                if(std::fabs(res_t) > tmx) tmx = std::fabs(res_t);
                n++;
            }
        }
        if(n > 0){
            std::cout << "      unbalanced acceleration left at iteration 0 (rms, non-dim):"
                      << std::endl
                      << "        radial      " << std::scientific << std::setprecision(3)
                      << sqrt(r2 / n) << "   (max " << rmx
                      << ", against " << sqrt(r2_0 / n) << " with no balance at all)"
                      << std::endl
                      << "        meridional  " << sqrt(t2 / n) << "   (max " << tmx
                      << ", against " << sqrt(t2_0 / n) << " with no balance at all)"
                      << std::fixed << std::endl;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for initBalancedState\n", elapsed.count() * 1e-9);
    std::cout << "      AGCM: initBalancedState end ......................." << std::endl;
}
