#ifndef _UTILS_
#define _UTILS_

#include <string>
#include <set>
#include <map>
#include <vector>

#include <iostream>
#include <fstream>
#include <limits>
#include <cstdint>
#include <cstring>
#include <cstdlib>   // getenv/atoi for coriolis_nontraditional()

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"

#define logger() \
if (false) ; \
else get_logger()

namespace AtomUtils{
    using namespace std;

    // ---- Which Coriolis approximation BOTH spheres use (ATOM_CORIOLIS_NONTRAD) ----
    //
    // The two used to disagree, silently, in a coupled model that exchanges momentum at the sea
    // surface. The full acceleration -2*Omega x u in (r, theta, phi) with theta = colatitude is
    //     F_r   = +2*Omega*sin(theta)*w
    //     F_the = +2*Omega*cos(theta)*w
    //     F_phi = -2*Omega*cos(theta)*v - 2*Omega*sin(theta)*u
    // RHS_Hyd_Turb carried all of it; RHS_Atm_Turb carried only the cos(theta) pair, dropping
    // F_r entirely and the sin(theta)*u part of F_phi. That truncation is the TRADITIONAL
    // approximation, and the atmosphere's own comment records why it got there: the +2*Omega*
    // sin(theta)*w term drove the eastward jets steadily upward under the lagging pressure
    // projection. It is also the textbook-correct way to truncate, because the two dropped terms
    // form one energetically consistent pair — dropping only one of them injects energy.
    //
    // DEFAULT 0 = traditional in BOTH spheres. That is the standard for hydrostatic, shallow
    // models (a 16 km atmosphere over a 200 m ocean is as shallow as it gets), it is what the
    // atmosphere already did, and it is the only choice consistent with the curvature terms this
    // model actually has — namely none: the uv/r, uw/r and cot(theta) metric terms are absent from
    // both advection operators, and the non-traditional Coriolis terms conserve energy only
    // together with them.
    //
    // ATOM_CORIOLIS_NONTRAD=1 restores the full non-traditional set in both spheres, which is the
    // hydrosphere's previous behaviour. This is an env knob and not a config entry on purpose:
    // atm and hyd are separate executables with separate XML files, and a duplicated setting is
    // exactly the kind of thing that drifts apart.
    inline bool coriolis_nontraditional(){
        static const bool v = [](){
            const char* e = getenv("ATOM_CORIOLIS_NONTRAD");
            return e ? (atoi(e) != 0) : false; }();
        return v;
    }

    // ---- The two missing halves of the spherical metric (both spheres, both default 0) ----
    //
    // ATOM_METRIC_CURVATURE — the curvature terms of the advection operator. transport_u/v/w are
    // the plain (u.grad) and the material derivative in spherical coordinates is not:
    //     D u_r  /Dt - (v^2 + w^2)/r
    //     D u_the/Dt + (u*v)/r - w^2*cot(theta)/r
    //     D u_phi/Dt + (u*w)/r + v*w*cot(theta)/r
    // Note this is NOT the shallow-atmosphere approximation, which drops uv/r and uw/r but KEEPS
    // the cot(theta) pair — and vw*cot(theta)/r is the meridian-convergence term that conserves
    // the zonal wind's angular momentum, i.e. the one that makes jets. ATJUP carries all of them
    // (RHS_Jup_Turb.cpp:405-411); these two spheres carry none.
    //
    // ATOM_METRIC_DIVERGENCE — the Poisson source. The solvers form
    //     du/dr + (1/r)dv/dthe + (1/(r sin)) dw/dphi
    // and the divergence in spherical coordinates also has +2u/r and +v*cot(theta)/r. Without
    // them the projection makes the field divergence-free with respect to an operator that is not
    // the divergence, so a real divergence survives it and the pressure is not the pressure that
    // would remove it.
    //
    // ⚠️ BOTH BELONG WITH A CORRECT METRIC RADIUS. Every term added here carries 1/r, and with the
    // default rad.z ~ 1.5 that r is ~400x too small, so the terms themselves come in ~400x too
    // strong. Pair them with ATM_METRIC_RADIUS=6370 (atmosphere).
    //
    // MEASURED at 20 iterations, so nobody has to guess how much that matters in practice:
    //   with the radius on   max relative change  u 3.6e-02, p_dyn 1.4e+00, everything else < 5e-04
    //   with the radius off  max relative change  u 4.8e-02, p_dyn 1.3e-01, nue 7.4e-02
    // So the wrong-radius case is percent-level rather than catastrophic — these terms are small
    // corrections, not dominant ones. p_dyn moves most in both cases because the divergence source
    // is what changed and p_dyn is small. That is a reason to pair them, not a reason to panic.
    //
    // The hydrosphere has no radius knob yet, so there these are for measurement only.
    inline bool metric_curvature(){
        static const bool v = [](){
            const char* e = getenv("ATOM_METRIC_CURVATURE");
            return e ? (atoi(e) != 0) : false; }();
        return v;
    }
    inline bool metric_divergence(){
        static const bool v = [](){
            const char* e = getenv("ATOM_METRIC_DIVERGENCE");
            return e ? (atoi(e) != 0) : false; }();
        return v;
    }

    struct HemisphereCoords{
        double lat, lon;
        string east_or_west, north_or_south;
    }; 

    HemisphereCoords convert_coords(double lon, double lat);

    std::ofstream& get_logger();

    std::tuple<double, int, int, int>
        max_diff(int i, int j, int k, const Array &a1, const Array &a2);

    // Bit-level non-finite check.  The Makefile uses -ffast-math, which implies
    // -ffinite-math-only; that lets the optimizer constant-fold std::isfinite(v)
    // to true, silently disabling every NaN/Inf guard built on it. IEEE-754 NaN
    // and ±Inf both have the 11 exponent bits all-set, so masking the raw bits
    // detects them regardless of -ffast-math. Use this everywhere you would
    // otherwise call std::isfinite / std::isnan / std::isinf.
    inline bool is_finite_safe(double v){
        std::uint64_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        return (bits & 0x7FF0000000000000ULL) != 0x7FF0000000000000ULL;
    }

    inline bool is_land(const Array& h, int i, int j, int k){
        return fabs(h.x[i][j][k] - 1) < std::numeric_limits<double>::epsilon();
    }
    inline bool is_air(const Array& h, int i, int j, int k){
        return !is_land(h,i,j,k);
    }
    inline bool is_water(const Array& h, int i, int j, int k){
        return !is_land(h,i,j,k);
    }
    inline bool is_ocean_surface(const Array& h, int i, int j, int k){
        return i==0 && !is_land(h, i, j, k);
    }
    inline bool is_land_surface(const Array& h, int i, int j, int k){
        return is_land(h, i, j, k) && !is_land(h, i+1, j, k);
    }
    inline bool is_east_coast(const Array& h, int j, int k){
        if(k == h.get_km()-1 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j, k+1);
    }
    inline bool is_west_coast(const Array& h, int j, int k){
        if( k == 0 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j, k-1);
    }
    inline bool is_north_coast(const Array& h, int j, int k){
        if( j == 0 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j-1, k);
    }
    inline bool is_south_coast(const Array& h, int j, int k){
        if(j == h.get_jm()-1 ) return false;//on grid boundary
        return is_land(h, 0, j, k) && !is_land(h, 0, j+1, k);
    }
    inline double parabola(double x){
        return x * x - 2.0 * x;
    }
    // x^2-2x mapped to [lower_bound, upper_bound]:
    // x=1 → lower_bound, x=0 or 2 → upper_bound
    inline double parabola_interp(double lower_bound, double upper_bound, double x){
        return upper_bound + (upper_bound - lower_bound) * (x * x - 2.0 * x);
    }
    inline double cosine_profile(double ratio){
    // ratio 0 -> Pol (Süd), ratio 1 -> Äquator, ratio 2 -> Pol (Nord)
    // Wir mappen ratio [0, 2] auf Cosinus-Argument [-PI/2, PI/2]
    // Ein einfacherer Weg für das gleiche Profil:
    // f(1) = 1 (Äquator), f(0) = 0 (Pol), f(2) = 0 (Pol)
    
//        return sin(M_PI_2 * ratio) * (2.0 - ratio); // Alternative Formel
    // ODER klassisch:
        return cos((ratio - 1.0) * M_PI_2) - 1.0; 
    }    
    inline double Agnesi(double a, double x){
        return pow(a, 3.0) // Versiera di Agnesi
            /(pow(a, 2.0) + pow(x, 2.0));
    }

    inline double exp_func(double T_K, const double co_1, const double co_2){
        return exp(co_1 * (T_K - 273.15) / (T_K - co_2));  // temperature in °K
    }


    int lon_2_index(double lon);
    int lat_2_index(double lat);
    int RunStart(string comment);
    int RunEnd(string comment, int Ma, int Run_start);

    double C_Dalton(int i, int j, int k, double coeff_Dalton, 
        double u_0, Array &u, Array &v, Array &w);
    double simpson(int n1, int n2, double dstep, Array_1D &value);
    double trapezoidal(int n1, int n2, double dstep, Array_1D &value);
    double rectangular(int n1, int n2, double dstep, Array_1D &value);
    double GetMean_3D(int jm, int km, Array &val_3D);
    double GetMean_2D(int jm, int km, Array_2D &val_2D);

    void read_IC(const string& fn, double** a, int jm, int km);
    void smooth_steps(int k, int im, int jm, Array &value);
    void smooth_tropopause(int jm, std::vector<double> &value);
    void load_map_from_file(const std::string& fn, std::map<float, float>& m);
    void CalculateNodeWeights(int jm, int km);
    void calculate_node_weights();
    //do the fft in all 3 directions
    void fft_gaussian_filter_3d(Array& data, int sigma);
    // do the fft in one direction
    // the valid directions are 'i', 'j' and 'k'
    // sigma is the Standard deviation for Gaussian kernel and controls how blurry the result will be
    //the result will replace the input array
    void fft_gaussian_filter_3d(Array& data, int sigma, char direction);
    void mirror_padding(double* data, size_t i_len, size_t p_len);
    void fft_gaussian_filter(double* _data, double* kernel, size_t len);
    //change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
    void move_data(double* data, int len);
    void move_data(std::vector<int>& data, int len);
    void move_data_to_new_arrays(int im, int jm, int km, double coeff, 
        std::vector<Array*>& arrays, std::vector<Array*>& new_arrays);
    void move_data_to_new_arrays( int jm, int km, double coeff, 
        std::vector<Array*>& arrays, 
        std::vector<Array*>& new_arrays, int i = 0);
    void use_presets(bool preset_1, bool preset_2, bool preset_3);





    // Shapiro (1-2-1) wiggle-damping filter for an Array field.
    //
    // Each enabled axis gets one filter pass per call to damp_wiggles():
    //   f_new[n] = f[n] + strength * 0.25 * (f[n-1] - 2*f[n] + f[n+1])
    //
    // Axis boundary conventions:
    //   k (longitude) — periodic wrap
    //   j (latitude)  — clamped (no-flux) at poles
    //   i (vertical)  — clamped at top/bottom
    //
    // i_surface: pointer to the 2D fluid-domain floor index array.
    //   Accepts either i_topography[j][k] (atmosphere: first cell above terrain)
    //   or i_bathymetry[j][k] (hydrosphere: first cell above the seafloor).
    //   In both cases cells with i < i_surface[j][k] are skipped (solid body).
    //   Pass nullptr to skip the mask and filter all cells.
    //
    // strength ∈ [0, 1]: 0 = no change, 1 = full Shapiro step.
    // passes: number of times the filter sequence is repeated.
    void damp_wiggles(Array& field,
                      const std::vector<std::vector<int>>* i_surface,
                      bool along_i, bool along_j, bool along_k,
                      double strength = 1.0,
                      int passes = 1);

    // Latitude-dependent zonal (φ) Shapiro filter — the classic lat-lon "polar
    // filter".  On a lat-lon grid the zonal cell width rm·sinθ·dφ shrinks toward
    // the poles, so the explicit zonal CFL number ~ |w|/(rm·sinθ·dφ) diverges and
    // short φ-wavelengths blow up at high latitude (root cause of the ATOM
    // post-155 polar instability).  This applies the same periodic 1-2-1 zonal
    // step as damp_wiggles along k only, but the number of passes grows toward the
    // poles so the effective zonal resolution is held to that at the reference
    // latitude.  Equatorward of lat_filter_start_deg nothing is touched.
    //
    //   passes(j) = clamp( round( (sinθ_ref/sinθ_j)² − 1 ), 0, max_passes )
    //   sinθ_j = sin(colat[j]) = cos(lat_j);  sinθ_ref = cos(lat_filter_start_deg)
    //
    // The square is not arbitrary: a 1-2-1 filter applied n times smooths over
    // ~√n grid cells, so holding the effective zonal resolution to the reference
    // latitude needs n ∝ (Δx_ref/Δx)² = (sinθ_ref/sinθ_j)².  A linear law left the
    // 60–72° band with 0 passes (round-down), which is where the residual 67°S
    // near-surface Antarctic-coast seed lived.  With start=40° the quadratic law
    // gives 67° ~3 passes and ramps to the cap by ~78°, while the storm-track band
    // equatorward of ~53° still gets 0 (large-scale Southern-Ocean eddies preserved).
    //
    // colat: pointer to the per-row colatitude in radians (Array_1D::z, i.e. the.z).
    //        Row count is taken from field.jm, so colat must have field.jm entries.
    // i_surface / boundary conventions (periodic k, solid-cell skip) match
    // damp_wiggles.  Pass nullptr for i_surface to filter all cells.
    // pass_gain multiplies the per-row pass count: passes[j] =
    // clamp(round(pass_gain*((sinθ_ref/sinθ_j)² − 1)), 0, max_passes). gain=1 is the
    // pure resolution-equalising law; gain>1 strengthens the high-mid-latitude band
    // (the quadratic is too flat at 50–60° to damp the high-lat near-surface mode there).
    void polar_zonal_filter(Array& field,
                            const double* colat,
                            const std::vector<std::vector<int>>* i_surface,
                            double lat_filter_start_deg = 40.0,
                            int    max_passes = 12,
                            double pass_gain = 1.0);

    // Orographic Shapiro filter — companion to polar_zonal_filter for the OTHER lat-lon
    // CFL hot-spot: steep topography at ANY latitude (Patagonian Andes 49°S/74°W, Atlas
    // 36°N/2°E, NZ Alps 46°S/170°E). At a steep cliff the near-surface horizontal velocity
    // develops a grid-scale (2Δ) oscillation that explicit advection amplifies into a CFL
    // blow-up; the polar filter keys on latitude and these sites are mid-latitude, so it
    // never reaches them. This applies `passes` 1-2-1 Shapiro sweeps along k (periodic)
    // AND j to `field`, but ONLY at air cells within `n_layers_above` of the surface in
    // columns whose surface index jumps by >= `steep_threshold` cells relative to a
    // horizontal neighbour. Solid neighbours use no-flux (current-cell substitution),
    // matching damp_wiggles / polar_zonal_filter. strength in [0,1] scales the coefficient
    // (0.25 at 1). All parameters are tunable. The cliff 2Δ mode lives in the lowest few
    // air cells, so n_layers_above is deliberately SMALL: with im~41 a deeper reach smears
    // ~12 of the 14 air levels above a high plateau (Tibet) and erases the orographic
    // deflection itself. The orography slope cap in init_topography removes most cliff
    // seeds at source, so this only needs a gentle near-surface touch.
    void orographic_shapiro_filter(Array& field,
                                   const std::vector<std::vector<int>>& i_surface,
                                   int    steep_threshold = 3,
                                   int    n_layers_above  = 3,
                                   int    passes          = 1,
                                   double strength        = 1.0);

    // Gentle global radial (i) Shapiro de-checkerboarding. The steep-orography velocity
    // blow-up is a 2Δ grid-scale checkerboard whose purest component is on the RADIAL axis
    // (wmode r=1.000), and it grows at near-surface cells that orographic_shapiro_filter's
    // steep-column mask never flags — so NO other filter (polar=k, orographic=j+k+steep)
    // reaches it. This applies `passes` 1-2-1 radial sweeps to ALL fluid cells: a single
    // 1-2-1 (coeff 0.25·strength) annihilates a level-to-level 2Δ oscillation in one pass
    // yet barely touches smooth vertical shear, so it is safe to run globally every iter.
    // No-flux (current-cell substitution) at the solid surface below; i=0/i=im-1 left to BCs.
    void radial_shapiro_filter(Array& field,
                               const std::vector<std::vector<int>>& i_surface,
                               int    passes   = 1,
                               double strength = 1.0);

    // Higher-order (4th-order Shapiro, n=2) radial de-checkerboarding. Response 1 - sin^4(kΔ/2):
    // still annihilates the 2Δ grid mode completely (factor 0 at k=π) but leaves resolved
    // vertical shear nearly intact (4Δ mode damped only 25% vs 50% for the 1-2-1) — so it can
    // run globally every iter WITHOUT spinning the jet/Hadley-Ferrel cells down the way the
    // 1-2-1 does. 5-point stencil; falls back to a gentle 1-2-1 adjacent to the surface/lid
    // where i±2 is unavailable. No-flux (current-cell substitution) into solid cells below.
    void radial_shapiro_filter_ho(Array& field,
                                  const std::vector<std::vector<int>>& i_surface,
                                  int    passes   = 1,
                                  double strength = 1.0);

    // Radial (i) de-checkerboarding for SCALAR fields, confined to orographic columns and
    // gated to vertical local extrema. Unlike radial_shapiro_filter_ho (global, every fluid
    // cell), this kills the single-cell 2Δz t/moisture checkerboard that seeds the NZ-Alps
    // blow-up (cold cell ~surface+7 over the steep slope-capped Alps;
    // [[project_nz_alps_t_checkerboard]]) WITHOUT touching the smooth bulk scalar field:
    //   (1) only steep/sloped columns (|Δ surface-index| ≥ steep_threshold to a neighbour),
    //   (2) only the first n_layers_above fluid cells of those columns,
    //   (3) only where the cell is a vertical local extremum ((vm−c)(vp−c) > 0) — the
    //       de-checkerboarding analogue of the minmod advection limiter, so a monotonic
    //       lapse rate / inversion is left bit-unchanged.
    // Deliberately UNLIKE the reverted 2026-06-08 *horizontal* scalar orographic_shapiro_filter
    // (which perturbed the global field and tipped the 62°N seam): radial axis + extremum gate
    // confine the action to actual checkerboard cells. No-flux (current-cell) into solid below.
    void orographic_radial_shapiro_filter(Array& field,
                                          const std::vector<std::vector<int>>& i_surface,
                                          int    steep_threshold = 2,
                                          int    n_layers_above  = 10,
                                          int    passes          = 2,
                                          double strength        = 1.0);

    // Soft steep-massif field smoothing — gentle horizontal (k then j) 1-2-1 Shapiro applied
    // to a field (u,v,w,p_dyn) ONLY over HIGH+STEEP contour columns (same mask as the init-
    // time terrain massif smoothing: i_surface ≥ high_thresh AND |Δi_surface|≥steep_thresh),
    // reaching n_layers_deep cells up the column. Cleans the locally unphysical field above
    // steep mountains (Himalaya) without touching plains/coasts/broad terrain/ocean or the
    // expected real flow. No-flux centre substitution at land/sub-surface neighbours.
    void steep_massif_field_smoothing(Array& field,
                                      const std::vector<std::vector<int>>& i_surface,
                                      int    high_thresh    = 6,
                                      int    steep_thresh   = 3,
                                      int    n_layers_deep  = 30,
                                      int    passes         = 1,
                                      double strength       = 1.0);

    // Near-surface coastal Rayleigh sponge — last-resort soft cap for the dry-seeded
    // velocity runaway at steep coasts ([[project-bc-atm-mc-extrap-overshoot]] /
    // [[project-orography-slope-cap]] Gulf-of-Alaska mode). The slope cap skips
    // ocean↔land cliffs (cannot raise the ocean side) and bcVelSurfSur is init-only;
    // without an in-loop friction term in RHS_Atm, coastal pressure-gradient force
    // drives a linear velocity growth at i=0–1 that saturates at the ±100 m/s BC cap.
    //
    // Mask: any (j,k) whose i_surface differs from a horizontal neighbour by >= 1 cell
    // (every coast plus any sloped column). Action: at the lowest `layers` air cells of
    // each flagged column, if |field| exceeds `threshold` (non-dim), multiply by
    // (1 − rate). Below threshold the cell is untouched, so normal coastal flow / sea
    // breeze / katabatic flow at < threshold·u_0 m/s is preserved. Above threshold the
    // value relaxes back exponentially with e-folding ≈ −1/ln(1−rate) iters.
    void coastal_velocity_sponge(Array& field,
                                 const std::vector<std::vector<int>>& i_surface,
                                 double threshold = 2.0,   // non-dim, =20 m/s at u_0=10
                                 double rate      = 0.05,
                                 int    layers    = 3);

    // Upper-troposphere Rayleigh sponge — relaxes zonal anomalies toward the zonal mean
    // in the topmost levels. The standard GCM remedy for the lat-lon grid CFL near the
    // poles AND the meridian seam at high altitude: as cell width rm·sinθ·dφ shrinks with
    // latitude AND the explicit timestep becomes marginal in the upper troposphere, short
    // zonal waves blow up. The polar zonal Shapiro filter handles the latitude axis at all
    // altitudes; this targets the COMPLEMENTARY high-altitude band where the seam (k=0/360)
    // mode reaches NaN first (observed 62°N / 10 km at iter 359).
    //
    // At each (i, j) compute the zonal mean of `field` (average over k), then relax each
    // cell toward that mean with strength α(i) = α_max · (i − i_start) / (im − 1 − i_start),
    // i.e. zero at i_start, full at the model top. Preserves the zonal mean exactly —
    // damps only zonal anomalies (where the CFL mode lives). Apply per RK4 step on u, v, w
    // after the other filters. α_max=0.2 gives ~5-iter e-folding of zonal anomalies at the
    // top, decaying to 0 at i_start.
    void upper_rayleigh_sponge(Array& field,
                               int    i_start   = 30,
                               double alpha_max = 0.2);

    // Extreme-peak removal filter for an Array field.
    //
    // Along each enabled axis the two direct neighbours of every cell define a
    // local background (their mean) and a variability scale (half-spread,
    // |f[n+1] - f[n-1]| / 2).  If the cell's deviation from the background
    // exceeds threshold * half_spread the cell is replaced by the background.
    // Smooth gradients are left intact; only isolated spikes are removed.
    //
    // Axis boundary conventions and i_surface masking are identical to
    // damp_wiggles(): k periodic, j/i clamped, solid cells skipped.
    //
    // threshold: sensitivity — lower values remove more; typical range 1–4.
    // passes:    number of times the scan is repeated per call.
    void remove_peaks(Array& field,
                      const std::vector<std::vector<int>>* i_surface,
                      bool along_i, bool along_j, bool along_k,
                      double threshold = 2.0,
                      int passes = 1);

    // 2-D overloads (Array_2D, axes j and k only).
    // Boundary conventions: k periodic, j clamped at poles.
    void damp_wiggles(Array_2D& field,
                      bool along_j, bool along_k,
                      double strength = 1.0,
                      int passes = 1);
    void remove_peaks(Array_2D& field,
                      bool along_j, bool along_k,
                      double threshold = 2.0,
                      int passes = 1);

    template<class T>
    void set_values(T* a, T value, int len){
        for(int i = 0 ; i < len; i++){
            a[i] = value;
        }
    }
    template<class T>
    void chessboard_grid(T** a, int j, int k, int jm, int km){
        T value_1 = 0.9, value_2 = 1.1;
        T* line_1 = new T[km];
        T* line_2 = new T[km];
        int i = 0;
        bool flip = true;
        while(i+k < km){
            if(flip){
                set_values(line_1+i, value_1, k);
                set_values(line_2+i, value_2, k);
            }else{
                set_values(line_1+i, value_2, k);
                set_values(line_2+i, value_1, k);
            }
            i+=k;
            flip = !flip;
        }
        if(i<km-1){
            if(flip){
                set_values(line_1+i, value_1, km-1-i);
                set_values(line_2+i, value_2, km-1-i);
            }else{
                set_values(line_1+i, value_2, km-1-i);
                set_values(line_2+i, value_1, km-1-i);
            }
        }
        int cnt=0;
        flip = true;
        for(i = 0; i<jm; i++){
            if(flip)
                memcpy(a[i], line_1, km*sizeof(T));
            else
                memcpy(a[i], line_2, km*sizeof(T));
            cnt++;
            if(cnt==j){
                flip = !flip;
                cnt = 0;
            }
        }
        delete[] line_1;
        delete[] line_2;
        return;
    }
}

#endif
