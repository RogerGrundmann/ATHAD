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
    constexpr double MAGNUS_A_WATER = 17.2694;
    constexpr double MAGNUS_B_WATER = 35.86;
    constexpr double MAGNUS_A_ICE = 21.8747;
    constexpr double MAGNUS_B_ICE = 7.66;

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
                const double E_wat = hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER);
                const double E_ice = hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE, MAGNUS_B_ICE);
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

                // Saturation specific humidity with safety fallback
                const double q_Satur = (p_u > E_sat) ? 
                                       (ep * E_sat / (p_u - E_sat)) : 
                                       (ep * FALLBACK_Q_FACTOR);

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
                if (p_u > E_sat_loc && E_sat_loc > 0.0) {
                    q_sat_loc = ep * E_sat_loc / (p_u - E_sat_loc);
                } else if (E_sat_loc <= 0.0) {
                    q_sat_loc = 0.0;
                } else {
                    q_sat_loc = ep * FALLBACK_Q_FACTOR;
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
                  
  cout << "tropopause_pole=" << tropopause_pole                                                                                                                                                            
       << " x_max=" << x_max                                                                                                                                                                               
       << " pole_index=" << round(tropopause_pole/L_atm) << endl;                                                                                                                                          



                                                                                                                                                                                         
    // Build symmetric cache of heights [m] and grid indices in one pass.                                                                                                                                
    std::vector<double> tropo_height_cache(jm);
    tropopause_layers = std::vector<double>(jm);                                                                                                                                                         
  
    for(int j = 0; j <= j_half; j++){                                                                                                                                                                    
        double x = x_max * (double)(j_half - j) / (double)j_half;
        double h = AtomUtils::Agnesi(tropopause_equator, x);                                                                                                                                             
        tropo_height_cache[j]       = h;
        tropo_height_cache[j_max-j] = h;                                                                                                                                                               
        tropopause_layers[j]        = round(h / L_atm);                                                                                                                                                 
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
    // Temperature Variables Declaration
    // ========================================================================
    double t_equat = 0.0; 
    double t_pole = 0.0; 
    double t_equat_add = 0.0; 
    double t_pole_add = 0.0; 
    double t_global_mean_exp = 0.0;
    double t_equat_curr = 0.0; 
    double t_pole_curr = 0.0; 
    double t_equat_prev = 0.0; 
    double t_pole_prev = 0.0;

    // ========================================================================
    // Step 1: Fix NASA Temperature Data Artifact at 180°E
    // ========================================================================
    if (is_first_time_slice()) {
        int k_half = (km - 1) / 2;

        // Interpolate bad data at dateline
        #pragma omp parallel for
        for (int j = 0; j < jm; j++) {
            temperature_NASA.y[j][k_half] = 
                0.5 * (temperature_NASA.y[j][k_half + 1] + 
                       temperature_NASA.y[j][k_half - 1]);              // [°C]
        }
    }

    // ========================================================================
    // Step 2: Initialize Temperature Field
    // ========================================================================
    const double inv_t0 = 1.0 / t_0;

    if (is_first_time_slice()) {
        // First time slice: Use NASA data directly
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < km; k++) {
            for (int j = 0; j < jm; j++) {
                t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0) * inv_t0;// Non-dimensional
            }
        }
    }

    // ========================================================================
    // Step 3: Apply EarthByte Reconstruction (if enabled)
    // ========================================================================
    if (!is_first_time_slice() && use_earthbyte_reconstruction) {
        #pragma omp parallel for collapse(2)
        for (int k = 0; k < km; k++) {
            for (int j = 0; j < jm; j++) {
                double val_nd = (t.x[0][j][k] + t_0) * inv_t0;
                t.x[0][j][k] = val_nd;                                  // Non-dimensional
                temp_reconst.y[j][k] = val_nd;                          // Store for later use
            }
        }
    }

    // ========================================================================
    // Step 4: Extract Temperature Values from Curves
    // ========================================================================
    t_equat_modern = get_temperatures_from_curve(0, m_equat_temperature_curve);
    t_pole_modern = get_temperatures_from_curve(0, m_pole_temperature_curve);
    t_global_mean_exp = get_temperatures_from_curve(*get_current_time(), 
                                                     m_global_temperature_curve);

    if (is_first_time_slice()) {
        t_global_mean_exp = get_temperatures_from_curve(0, m_global_temperature_curve);
        t_global_mean = GetMean_2D(jm, km, temperature_NASA);
    }

    t_global_mean = get_temperatures_from_curve(*get_current_time(), 
                                                 m_global_temperature_curve);

    // ========================================================================
    // Step 5: Calculate Temperature Increments (for non-first time slices)
    // ========================================================================
    // Current-Ma equatorial/polar temperatures from the Scotese curves. These
    // ALONE define the parabolic paleo profile applied in Step 7 — there is no
    // dependence on any foregoing Ma. Computed for every paleo slice, whether or
    // not a preceding slice exists, so a single-Ma run (time_start == time_end)
    // works without prior-slice reconstruction.
    if (*get_current_time() > 0) {
        t_equat = get_temperatures_from_curve(*get_current_time(), m_equat_temperature_curve);
        t_pole  = get_temperatures_from_curve(*get_current_time(), m_pole_temperature_curve);
    }

    // The inter-slice increments below feed ONLY the optional EarthByte
    // reconstruction correction and need a preceding slice — get_previous_time()
    // throws on the first slice — so keep them guarded by !is_first_time_slice().
    if (!is_first_time_slice()) {
        // Temperature changes between time steps
        t_equat_add = get_temperatures_from_curve(*get_current_time(), m_equat_temperature_curve)
                    - get_temperatures_from_curve(*get_previous_time(), m_equat_temperature_curve);

        t_pole_add = get_temperatures_from_curve(*get_current_time(), m_pole_temperature_curve)
                   - get_temperatures_from_curve(*get_previous_time(), m_pole_temperature_curve);

        // Current and previous values
        t_equat_curr = get_temperatures_from_curve(*get_current_time(), m_equat_temperature_curve);
        t_equat_prev = get_temperatures_from_curve(*get_previous_time(), m_equat_temperature_curve);
        t_pole_curr = get_temperatures_from_curve(*get_current_time(), m_pole_temperature_curve);
        t_pole_prev = get_temperatures_from_curve(*get_previous_time(), m_pole_temperature_curve);
    }

    // ========================================================================
    // Step 6: Print Diagnostics
    // ========================================================================
    std::cout.precision(3);
    std::cout << "\n       Time slice of Paleo-AGCM: ...................... Ma = " << Ma << " million years\n";
    std::cout << "\n       Equatorial temperature increase: ................ t_equat_add      = " << t_equat_add << " °C";
    std::cout << "\n       Polar temperature increase: ..................... t_pole_add       = " << t_pole_add << " °C";
    std::cout << "\n       Equatorial temperature at paleo times: .......... t_equat_paleo    = " << t_equat << " °C";
    std::cout << "\n       Polar temperature at paleo times: ............... t_pole_paleo     = " << t_pole << " °C";
    std::cout << "\n       Mean temperature at paleo times: ................ t_global_mean    = " << t_global_mean << " °C";
    std::cout << "\n       Expected mean temperature at paleo times: ....... t_global_mean_exp= " << t_global_mean_exp << " °C";
    std::cout << "\n       Equatorial temperature at modern times: ......... t_modern_equat   = " << t_equat_modern << " °C";
    std::cout << "\n       Polar temperature at modern times: .............. t_modern_pole    = " << t_pole_modern << " °C\n\n";

    // ========================================================================
    // Step 7: Apply Latitudinal Temperature Distribution
    // ========================================================================
    const double d_j_half = 0.5 * (jm - 1);
    const double t_0_inv = 1.0 / t_0;

    // Convert to non-dimensional
    t_equat = (t_equat + t_0) * t_0_inv;
    t_pole = (t_pole + t_0) * t_0_inv;
    t_equat_add = (t_equat_add + t_0) * t_0_inv;
    t_pole_add = (t_pole_add + t_0) * t_0_inv;

    // Effective temperature gradients
    const double delta_equat_nd = (t_equat_curr - t_equat_prev) / t_0;
    const double delta_pole_nd = (t_pole_curr - t_pole_prev) / t_0;
    const double delta_t_eff = delta_pole_nd - delta_equat_nd;
    const double t_eff = t_pole - t_equat;

    // Modern slice (Ma == 0) retains the observed NASA field assigned in Step 2;
    // paleo slices (Ma > 0) get the Scotese pole→equator parabola.
    const bool modern = (*get_current_time() == 0);

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            double ratio = (double)j / d_j_half;

            if (!use_earthbyte_reconstruction) {
                // Standard parabolic pole-to-pole distribution, built solely from
                // the current Ma's Scotese equator/pole temperatures — no foregoing
                // Ma required. Ma == 0 keeps the NASA field set in Step 2.
                if (!modern) {
                    t.x[0][j][k] = t_eff * AtomUtils::parabola(ratio) + t_pole;
                }
            } else {
                // EarthByte reconstruction active
                if (*get_current_time() == 0) {
                    // Initial state from NASA data
                    t.x[0][j][k] = (temperature_NASA.y[j][k] + t_0) * t_0_inv;
                } else {
                    // Apply correction to maintain latitudinal gradient evolution
                    double correction_nd = delta_t_eff * AtomUtils::parabola(ratio) + delta_pole_nd;
                    t.x[0][j][k] = temp_reconst.y[j][k] + correction_nd;
                }
            }
        }
    }

    // ========================================================================
    // Step 8: Vertical Temperature Profile & Potential Temperature
    // ========================================================================
    const double beta = BETA_COSMO;
    const double R_W_R_A = R_WaterVapour / R_Air;

    #pragma omp parallel for collapse(2)
    for (int k = 0; k < km; k++) {
        for (int j = 0; j < jm; j++) {
            int i_mount = i_topography[j][k];
            double t_u_init = t.x[0][j][k] * t_0;                       // [K]

            // OPTION B: broad land-sea thermal contrast (stacks on top of Option A).
            // The high-heat-capacity ocean stays near the zonal parabola; land equilibrates
            // closer to radiative equilibrium -> WARMER than the zonal reference in low
            // latitudes (subtropical/tropical continents -> thermal lows / monsoons) and
            // COLDER at high latitudes (continental interiors). cos(2*lat) gives +amp at the
            // equator, 0 near 45 deg, -amp at the poles. Applied to the SEA-LEVEL reference,
            // so Option A's elevation cooling still stacks on top. Paleo land only; the modern
            // NASA field and all ocean cells are untouched.
            constexpr double LANDSEA_AMP = 8.0;                        // [K] land-sea contrast amplitude (prototype, tunable)
            if (*get_current_time() != 0 && is_land(h, 0, j, k)) {
                const double lat_rad = (90.0 - (double)j * 180.0 / (double)(jm - 1)) * M_PI / 180.0;
                t_u_init += LANDSEA_AMP * cos(2.0 * lat_rad);
            }

            // Surface-temperature anchor.
            // Modern (Ma==0): the NASA field is observed AT the terrain top, so project it
            // dry-adiabatically to sea level; the COSMO column build below then reproduces it
            // exactly at i_mount (legacy behaviour, unchanged).
            // Paleo (Ma>0, OPTION A): the Scotese parabola is the SEA-LEVEL latitudinal
            // reference, NOT the mountain-top value. Anchor it at i=0 WITHOUT projecting, so
            // the COSMO vertical profile cools each column by its own DEM elevation:
            //     t.x[i_mount] = sqrt(T_sl^2 - 2*beta*g*h_mount/R)      (~5.5 K/km)
            // -> cold plateaus / warm lowlands -> zonal thermal structure that restores the
            // topographically-anchored baroclinic eddies the zonally-constant profile killed.
            // Ocean columns (i_mount=0) carry no elevation, so they stay at the parabola.
            if (*get_current_time() == 0) {                             // modern: project mountain-top NASA value up to sea level
                auto [t_pot, p_stat_0] = project_to_sea_level(
                    t_u_init,
                    get_layer_height(i_mount),
                    beta, R_Air, r_air, g
                );
                t.x[0][j][k] = t_pot;
                p_stat.x[0][j][k] = p_stat_0;
            } else {                                                    // paleo OPTION A: parabola IS the sea-level reference (no projection)
                p_stat.x[0][j][k] = 1e-2 * (r_air * R_Air * t_u_init);
                t.x[0][j][k]      = t_u_init;
            }

            temp_pot.y[j][k]     = t.x[0][j][k];
            temp_reconst.y[j][k] = t.x[0][j][k];

            // ================================================================
            // Surface Humidity Calculation (Magnus Formula)
            // ================================================================
            const double t_u_0 = t_u_init;
            const double p_u_0 = p_stat.x[0][j][k];
            
            // Magnus coefficients depend on phase (water vs ice)
            const double a_loc = (t_u_0 >= t_0) ? MAGNUS_A_WATER : MAGNUS_A_ICE;
            const double b_loc = (t_u_0 >= t_0) ? MAGNUS_B_WATER : MAGNUS_B_ICE;

            const double E_Satur = hp * exp(a_loc * (t_u_0 - t_0) / (t_u_0 - b_loc));
            const double e_curr  = c.x[0][j][k] * p_u_0 / (c.x[0][j][k] + ep);

            relative_humidity.y[j][k] = std::clamp(
                (e_curr / E_Satur) * 100.0,
                MIN_HUMIDITY,
                MAX_HUMIDITY
            );

            // ================================================================
            // Vertical Profile: Temperature, Pressure, Density
            // ================================================================
            const double t_safe = std::max(MIN_SAFE_TEMP, t.x[0][j][k]);
            const double tu_be     = t_safe / beta;
            const double c_inv_tu2 = (2.0 * beta * g) / (R_Air * t_safe * t_safe);
            const double inv_R_Air = 1.0 / R_Air;
            const double p_basis   = p_stat.x[0][j][k];

            // Removed inner #pragma omp simd to avoid nested parallelization issues
            for (int i = 0; i < im; i++) {
                const double height = get_layer_height(i);
                const double s_i    = sqrt(std::max(0.0, 1.0 - height * c_inv_tu2));

                const double t_curr = t_safe * s_i;  // [K]
                const double p_val  = p_basis * exp(-tu_be * (1.0 - s_i));

                // Smooth lower bound: approaches t_00 asymptotically rather than clamping hard
                const double delta = t_curr - t_00;
                const double sharpness = 5.0;                           // larger = closer to hard clamp

                t.x[i][j][k] = t_00 + delta / (1.0 + std::exp(-sharpness * delta));
                t.x[i][j][k] = std::max(t_00, t_curr);                  // Ensure minimum temperature
                p_stat.x[i][j][k] = p_val;

                // Density calculations
                const double rho_base   = (p_val * 100.0 * inv_R_Air);
                const double inv_t_curr = 1.0 / std::max(MIN_SAFE_TEMP, t_curr);

                const double virtual_mult = 1.0 + (R_W_R_A - 1.0) * c.x[i][j][k];
                const double total_water_factor = 
                    std::max(MIN_WATER_FACTOR, 1.0 - (cloud.x[i][j][k] + ice.x[i][j][k]));

                const double mask = is_land(h, i, j, k) ? r_air : 1.0;

                r_dry.x[i][j][k]   = (rho_base * inv_t_curr) * mask;
                r_humid.x[i][j][k] = (rho_base / (t_curr * virtual_mult * total_water_factor)) * mask;
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
*
*/
void cAtmosphereModel::load_global_temperature_curve(){
    load_map_from_file(temperature_global_file, m_global_temperature_curve);
/*
    cout << "   m_global_temperature_curve" << endl;
    for(const auto &printout : m_global_temperature_curve){
        cout << printout.first << " ..... " << printout.second << '\n';
    }
*/
}
/*
*
*/
void cAtmosphereModel::load_equat_temperature_curve(){
    load_map_from_file(temperature_equat_file, m_equat_temperature_curve);
/*
    cout << "   m_equat_temperature_curve" << endl;
    for(const auto &printout : m_equat_temperature_curve){
        cout << printout.first << " ..... " << printout.second << '\n';
    }
*/
}
/*
*
*/
void cAtmosphereModel::load_pole_temperature_curve(){
    load_map_from_file(temperature_pole_file, m_pole_temperature_curve);
/*
    cout << "   m_pole_temperature_curve" << endl;
    for(const auto &printout : m_pole_temperature_curve){
        cout << printout.first << " ..... " << printout.second << '\n';
    }
*/
}
/*
*
*/
float cAtmosphereModel::get_temperatures_from_curve(float time, 
    std::map<float, float>& m) const{
    if(time < m.begin()->first 
        || time > (--m.end())->first){
        std::cout << "Input time out of range: " << time << std::endl;    
        return NAN;
    }
    if(m.size() < 2){
        std::cout << "No enough data in map m" << std::endl;
        return NAN;
    }
    map<float, float>::const_iterator upper = m.begin(), 
        bottom = ++m.begin(); 
    for(map<float, float>::const_iterator it = m.begin();
            it != m.end(); ++it){
        if(time < it->first){
            bottom = it;
            break;
        }else{
            upper = it;
        }
    }
/*
    std::cout << "   get_temperatures_from_curve" << std::endl;
    std::cout << "   Ma ->   " << upper->first << " " 
        << bottom->first << std::endl;
    std::cout << "   temp-range ->   "<< upper->second 
        << " " << bottom->second << std::endl;
    std::cout << "   temp-interpolation ->   " 
        << upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second) << std::endl << std::endl;
*/
    return upper->second + (time - upper->first) 
       /(bottom->first - upper->first) 
        * (bottom->second - upper->second);
}
/*
*
*/
void cAtmosphereModel::LandOceanFraction(){
// calculation of the ratio ocean to land, also addition and subtraction of CO2 of land, ocean and vegetation

    cout << endl << endl << endl << "      AGCM: LandOceanFraction" << endl;

    int h_point_max = (jm-1) * (km-1);
    int h_land = 0;

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            if(is_land(h, 0, j, k))  h_land = h_land + h.x[0][j][k];
        }
    }

    int h_ocean = h_point_max - h_land;
    double ocean_land = (double)h_ocean/(double)h_land;

    cout.precision(3);
    cout << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      total number of points at constant height " << " = " 
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') 
        << h_point_max << endl << setiosflags(ios::left) << setw(50) 
        << setfill('.') << "      number of points on the ocean surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << h_ocean << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      number of points on the land surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << h_land << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      ocean/land ratio " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << ocean_land 
        << endl << endl;
    cout << setiosflags(ios::left) << setw(50) << setfill('.') 
        << "      addition of CO2 by ocean surface " << " = " 
        << resetiosflags(ios::left) << setw(7) << fixed << setfill(' ') 
        << co2_ocean << endl << setiosflags(ios::left) << setw(50) 
        << setfill('.') << "      addition of CO2 by land surface " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << co2_land << endl << setiosflags(ios::left) 
        << setw(50) << setfill('.') << "      subtraction of CO2 by vegetation " 
        << " = " << resetiosflags(ios::left) << setw(7) << fixed 
        << setfill(' ') << co2_vegetation << endl << setiosflags(ios::left) 
        << setw(50) << "      valid for one single point on the surface"<< endl << endl;
    cout << endl;


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

            const double RH_init = is_land(h, i_mount, j, k) ? 0.60 : 0.75;
            for (int i = 0; i < im; i++) {
                double t_u = t.x[i][j][k] * t_0;
                double p_u = p_stat.x[i][j][k];

                const double E_sat = (t_u >= t_0)
                    ? hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER)
                    : hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE,   MAGNUS_B_ICE);

                const double q_sat = (p_u > E_sat) ? ep * E_sat / (p_u - E_sat)
                                                    : ep * FALLBACK_Q_FACTOR;

                c.x[i][j][k]     = (i >= i_mount) ? RH_init * q_sat : 0.0;
//                c.x[i][j][k]     = 1.5 * c.x[i][j][k];                  // a very big cloud at 2 km and tends to reach the ground in higher latitudes
                c.x[i][j][k]     = 1.25 * c.x[i][j][k];                 // gives a nice cloud around 1 km height
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
    // Parabola H_crit(p): roots at p=0 and p=p_crit, minimum Hu_cr_mid at p=p_mid.
    // H_crit = Hu_cr_max - Hu_curv * x * (1 - x),  x = p / p_crit
    const double p_crit     = 1000.0;
    const double p_mid      = 550.0;
    const double x_mid      = p_mid / p_crit;                           // 0.55
    const double Hu_curv    = Hu_diff / (x_mid * (1.0 - x_mid));        // ~0.808
    const double inv_p_crit = 1.0 / p_crit;

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

                const double E_sat = (t_u >= t_0)
                    ? hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER)
                    : hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE,   MAGNUS_B_ICE);
                const double q_sat = ep * E_sat / (p_u - E_sat);

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
            for (int i = 0; i < im; i++) {
                const double t_u = t.x[i][j][k] * t_0;
                const double p_u = p_stat.x[i][j][k];

                const double E_sat = (t_u >= t_0)
                    ? hp * AtomUtils::exp_func(t_u, MAGNUS_A_WATER, MAGNUS_B_WATER)
                    : hp * AtomUtils::exp_func(t_u, MAGNUS_A_ICE,   MAGNUS_B_ICE);
                const double q_sat = ep * E_sat / (p_u - E_sat);

                const double x_norm = p_u * inv_p_crit;
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
