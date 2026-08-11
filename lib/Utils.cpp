#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <algorithm>

#include <Utils.h>
#include "cAtmosphereModel.h"

using namespace AtomUtils;
using namespace std;

std::vector<std::vector<double> > m_node_weights;

// Cosine-of-latitude area weights, built once on first use.
//
// These used to be built lazily inside GetMean_2D/GetMean_3D by the sequence
// "if(size != jm) CalculateNodeWeights()", which clears and re-push_backs the global
// vector. printDataAtm() calls those means from ~9 concurrent #pragma omp sections, so
// on the FIRST call every thread saw the size mismatch and raced into the same clear()
// + push_back() — reallocating one vector from several threads at once and corrupting
// the heap ("double free or corruption"). It was latent in ATOM_Precipitation only
// because an earlier single-threaded GetMean_2D(temperature_NASA) call in
// initTemperatureData happened to populate the weights first; ATHAD has no NASA field,
// so printDataAtm became the first caller and the race went live.
//
// A function-local static is initialised exactly once, and concurrent callers block
// until that completes (C++11 [stmt.dcl]/4), which is precisely the guarantee needed.
static const std::vector<std::vector<double> >& node_weights(int jm, int km){
    static const std::vector<std::vector<double> > w = [jm, km]{
        std::vector<std::vector<double> > v;
        v.reserve(jm);
        for(int i = 0; i < jm; i++){
            const double weight = (i <= 90) ? cos((90 - i) * M_PI / 180.0)
                                            : cos((i - 90) * M_PI / 180.0);
            v.emplace_back(km, weight);
        }
        return v;
    }();
    return w;
}

std::ofstream& AtomUtils::get_logger(){
    static std::ofstream o("atom_log.txt", std::ofstream::out);
    return o;
}
/*
*
*/
HemisphereCoords AtomUtils::convert_coords(double lon, double lat){
    HemisphereCoords ret;
    if ( lat > 90 ){
        ret.lat = lat - 90;
        ret.north_or_south = "°S";
    }else{
        ret.lat = 90 - lat;
        ret.north_or_south = "°N";
    }
    if ( lon > 180 ){
        ret.lon = 360 - lon;
        ret.east_or_west = "°W";
    }else{
        ret.lon = lon;
        ret.east_or_west = "°E";
    }
    return ret;
}
/*
* [-180, 180] to [0, 360]
*/
int AtomUtils::lon_2_index(double lon){
    if(lon<0)
        return int(round(360+lon));
    else
        return int(round(lon));
}
/*
* [-90, 90] to [0, 180]
*/
int AtomUtils::lat_2_index(double lat){
    return int(round(90-lat));
}
/*
*
*/
//change data coordinate system from -180° _ 0° _ +180° to 0°- 360°
void AtomUtils::move_data(double* data, int len){
    for(int i=0; i<len/2; i++){
        std::iter_swap(data+i, data+len/2+i);
    }
    data[len-1] = data[0];
    return;
}
/*
*
*/
void AtomUtils::move_data(std::vector<int>& data, int len){   
    for(int i=0; i<len/2; i++){
        std::iter_swap(data.begin()+i, data.begin()+len/2+i);
    }
    data[len-1] = data[0];
    return;
}
/*
*
*/
//the size of arrays must be the same with the size of new_arrays
//the values in the new_arrays will be assigned to coeff * old_values
void AtomUtils::move_data_to_new_arrays(int im, int jm, int km, double coeff, 
    std::vector<Array*>& arrays, std::vector<Array*>& new_arrays){
    assert(arrays.size() == new_arrays.size());
    for ( int i = 0; i < im; i++ ){
        move_data_to_new_arrays(jm, km, coeff, arrays, new_arrays, i);
    }
    return;
}
/*
*
*/
void AtomUtils::move_data_to_new_arrays(int jm, int km, double coeff, 
    std::vector<Array*>& arrays, std::vector<Array*>& new_arrays, int i){
    assert(arrays.size() == new_arrays.size());
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++ ){
            for(std::size_t n=0; n<arrays.size(); n++){
                new_arrays[n]->x[ i ][ j ][ k ] = coeff * arrays[n]->x[ i ][ j ][ k ];
            }
        }
    }
    return;
}
/*
*
*/
std::tuple<double, int, int, int>
AtomUtils::max_diff(int im, int jm, int km, const Array &a1, const Array &a2){
    std::tuple<double, int, int, int> ret;
    for ( int i = 0; i < im; i++ ){
        for ( int j = 0; j < jm; j++ ){
            for ( int k = 0; k < km; k++ ){
                double tmp = fabs ( a1.x[ i ][ j ][ k ] - a2.x[ i ][ j ][ k ] );
                if(tmp > std::get<0>(ret)){
                    std::get<0>(ret) = tmp;
                    std::get<1>(ret) = i;
                    std::get<2>(ret) = j;
                    std::get<3>(ret) = k;
                }
            }
        }
    }
    return ret;
}
/*
*
*/
double AtomUtils::C_Dalton(int i, int j, int k, double coeff_Dalton, 
    double u_0, Array &u, Array &v, Array &w){
    // variation of the evaporation coefficient in Dalton's evaporation law, parabola
    // air velocity measured 2m above sea level
    // for v_max = 10 m/s, but C_Dalton is function of v, should be included
    // Geiger ( 1961 ) by > Zmarsly, Kuttler, Pethe in mm/( h * hPa ), p. 133
    // for the modern world the global precipitation is 10% higher than evaporation

    double v_max = 10.0;                                                // in [m/s] ..... Geiger ( 1961 ) by Zmarsly, Kuttler, Pethe in m/s, p. 133
    double c_Dalton_max = 0.053;                                        // in [mm/(h*hPa)]
    double vel_magnitude = sqrt((u.x[i][j][k] * u.x[i][j][k] 
                               + v.x[i][j][k] * v.x[i][j][k] 
                               + w.x[i][j][k] * w.x[i][j][k]) / 3.0) * u_0;// in [m/s]
    return coeff_Dalton 
        * sqrt(c_Dalton_max * c_Dalton_max * vel_magnitude/v_max);      // result in [mm/(h*hPa)]
}
/*
*
*/
void AtomUtils::read_IC(const string& fn, double** a, int jm, int km){

    cout << endl << "      read_IC" << endl;

    cout << "      " << fn << endl;

    ifstream ifs(fn);

    if(!ifs.is_open()){
        cerr << "ERROR: unable to open " << fn << "\n";
        abort();
    }

    double lat, lon, d;

    for(int k=0; k < km && !ifs.eof(); k++){
        for(int j=0; j < jm; j++){
            ifs >> lat >> lon >> d;
            a[j][k] = d;
        }
    }

    cout << "      read_IC ended" << endl;
    return;
}
/*
*
*/
void AtomUtils::smooth_steps(int k, int im, int jm, Array &value){
    double coeff = 0.5;
    std::vector<std::vector<double> > inter(im, std::vector<double>(jm, 0));
    for(int i = 0; i < im ; i++){
        for(int j = 0; j < jm; j++){
            inter[i][j] = value.x[i][j][k];
        }
    }
    for(int i = 1; i < im-1 ; i++){
        for(int j = 1; j < jm-1; j++){
            value.x[i][j][k] = (coeff * inter[i][j] 
                + (1.0 - coeff) * (inter[i+1][j] 
                + inter[i-1][j] + inter[i][j+1] 
                + inter[i][j-1])/4.0);
        } 
    }
    return;
}
/*
*
*/
void AtomUtils::smooth_tropopause(int jm, std::vector<double> &value){
    double coeff = 0.5;
    std::vector<double>inter(jm, 0);
    for(int j = 0; j < jm; j++){
        inter[j] = value[j];
    }
    for(int j = 2; j < jm-2; j++){
        value[j] = (coeff * inter[j] + (1.0 - coeff) 
            * (inter[j+1] + inter[j-1] + inter[j+2] + inter[j-2])/4.0);
    } 
    return;
}
/*
*
*/
double AtomUtils::simpson(int n1, int n2, double dstep, Array_1D &value){
        double sum_even=0, sum_odd=0;
        if(n2 % 2 == 0){
            for(int i = n1+1; i < n2; i+=2){sum_odd += 4 * value.z[i];}
            for(int i = n1+2; i < n2; i+=2){sum_even += 2 * value.z[i];}
        }else cout << "       n2    must be an even number to use the \
            Simpson integration method" << endl;
    return dstep/3 * (value.z[n1] + sum_odd + sum_even + value.z[n2]); // Simpson Rule integration
}
/*
*
*/
double AtomUtils::trapezoidal(int n1, int n2, double dstep, Array_1D &value){
        double sum = 0;
        for(int i = n1+1; i < n2; i++){sum += 2 * value.z[i];}
        return dstep/2 * (value.z[n1] + sum + value.z[n2]);  // Trapezoidal Rule integration
}
/*
*
*/
double AtomUtils::rectangular(int n1, int n2, double dstep, Array_1D &value){
        double sum = 0;
        for(int i = n1; i <= n2; i++){sum += value.z[i];}
        return dstep * sum;  // Rectangular Rule integration
}
/*
*
*/
void AtomUtils::load_map_from_file(const std::string& fn, 
                                   std::map<float, float>& m) {
    std::ifstream f(fn);
    if (!f.is_open()) {
        std::cerr << "error occurred while opening file: " << fn << std::endl;
        return;
    }

    std::string line;
    while (std::getline(f, line)) {
//        if (line.empty() || line[0] == '#') continue;
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        float key = 0, val = 0;
        std::stringstream ss(line);
        
        if (ss >> key >> val) {
            m[key] = val;
        }
    }
}
/*
*
*/
double AtomUtils::GetMean_3D(int jm, int km, Array &val_3D){
    const std::vector<std::vector<double> >& nw = node_weights(jm, km);
    double ret=0.0, weight=0.0;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            ret += val_3D.x[0][j][k] * nw[j][k];
            weight += nw[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
double AtomUtils::GetMean_2D(int jm, int km, Array_2D &val_2D){
    const std::vector<std::vector<double> >& nw = node_weights(jm, km);
    double ret=0.0, weight=0.0;
    for(int j=0; j<jm; j++){
        for(int k=0; k<km; k++){
            ret += val_2D.y[j][k] * nw[j][k];
            weight += nw[j][k];
        }
    }
    return ret/weight;
}
/*
*
*/
void AtomUtils::CalculateNodeWeights(int jm, int km){
    // Kept as the public entry point, but it no longer mutates shared state: it just
    // forces the one-time, thread-safe build in node_weights(). Callers may still invoke
    // it eagerly before a parallel region; they no longer have to.
    (void)node_weights(jm, km);
    return;
}
/*
*
*/
int AtomUtils::RunStart(string comment){
    string at = "AGCM";
    string hy = "OGCM";
    string comment_1 = "";
    string comment_2 = "";
    if(comment.compare(at) == 0){
        comment_1 = " ... AGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
             Atmosphere General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    if(comment.compare(hy) == 0){
        comment_1 = " ... OGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
            Ocean General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    std::time_t Run_start;
    struct tm * timeinfo_begin;
    std::time(&Run_start);
    timeinfo_begin = std::localtime(&Run_start);
    std::cout << std::endl << std::endl;
    std::cout << comment_2 
        << std::endl << std::endl;
    std::cout << comment_1 
        << std::asctime(timeinfo_begin);
    return Run_start;
}
/*
*
*/
int AtomUtils::RunEnd(string comment, int Ma, int Run_start){
    string at = "AGCM";
    string hy = "OGCM";
    string comment_1 = "";
    string comment_2 = "";
        if(comment.compare(at) == 0){
        comment_1 = " ... AGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
             Atmosphere General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    if(comment.compare(hy) == 0){
        comment_1 = " ... OGCM: time and date at run time begin:   ";
        comment_2 = "    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \
            Ocean General Circulation Model \
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n";
    }
    std::time_t Run_end;
    std::time(&Run_end);
    struct tm * timeinfo_end;
    timeinfo_end = std::localtime(&Run_end);
    std::cout << std::endl << std::endl;
    std::cout << comment_2 
        << std::endl << std::endl;
    std::cout << comment_1 
        << std::asctime(timeinfo_end);
    int Run_total = Run_end - Run_start;
    int Run_total_minutes = Run_total/60;
    int Run_total_hours = Run_total_minutes/60;
    std::cout << std::endl << " ... computed time slice:"
        << "  Ma = " << Ma << std::endl << std::endl
        << " ... computer time needed:" << std::endl << std::endl
        << setw(20) << setfill(' ') << Run_total 
        << " seconds" << std::endl
        << " ... compares to:" << std::endl << std::endl
        << setw(20) << setfill(' ') << Run_total_hours 
        << " hours" << std::endl
        << setw(20) << setfill(' ') << Run_total_minutes 
        << " minutes" << std::endl
        << setw(20) << setfill(' ') << Run_total%60 
        << " seconds" << std::endl << std::endl
        << std::endl;
    return 0;
}
/*
*
*/
void AtomUtils::use_presets(bool preset_1, bool preset_2, bool preset_3){
    if(preset_1)
        cout << endl << "      use_earthbyte_reconstruction is in use" << endl;
    else
        cout << endl << "      use_earthbyte_reconstruction is NOT in use" << endl;

    if(preset_2)
        cout << "      use_NASA_temperature is in use" << endl;
    else
        cout << "      use_NASA_temperature is NOT in use" << endl;

    if(preset_3)
        cout << "      use_NASA_velocity is in use" << endl;
    else
        cout << "      use_NASA_velocity is NOT in use" << endl;
    return;
}

// ============================================================================
// Shapiro (1-2-1) wiggle-damping filter
// ============================================================================
void AtomUtils::damp_wiggles(Array& field,
                              const std::vector<std::vector<int>>* i_surface,
                              bool along_i, bool along_j, bool along_k,
                              double strength,
                              int passes)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    // Coefficient for the Shapiro update: f += coeff*(f_{-1} - 2f + f_{+1})
    const double coeff = strength * 0.25;

    // Helper: true when cell (i,j,k) is in the fluid domain.
    // Works for both i_topography (atmosphere) and i_bathymetry (hydrosphere):
    // in both cases i_surface[j][k] is the first fluid cell; solid cells have i < i_surface[j][k].
    auto in_fluid = [&](int i, int j, int k) -> bool {
        if (!i_surface) return true;
        // std::max guards against i_surface[j][k] < 0 (bathymetry deeper than L_hyd)
        return i >= std::max((*i_surface)[j][k], 0);
    };

    // Scratch buffer — snapshot before each axis pass to avoid in-place corruption.
    std::vector<double> tmp(im * jm * km);

    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        const int km1 = (k - 1 + km) % km;
                        const int kp1 = (k + 1)      % km;
                        // no-flux at solid neighbours: substitute current cell value
                        const double v_km1 = in_fluid(i, j, km1) ? tmp[idx(i,j,km1)] : tmp[idx(i,j,k)];
                        const double v_kp1 = in_fluid(i, j, kp1) ? tmp[idx(i,j,kp1)] : tmp[idx(i,j,k)];
                        field.x[i][j][k] = tmp[idx(i,j,k)]
                            + coeff * (v_km1 - 2.0 * tmp[idx(i,j,k)] + v_kp1);
                    }
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    const int jm1 = std::max(j - 1, 0);
                    const int jp1 = std::min(j + 1, jm - 1);
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        // no-flux at solid neighbours: substitute current cell value
                        const double v_jm1 = in_fluid(i, jm1, k) ? tmp[idx(i,jm1,k)] : tmp[idx(i,j,k)];
                        const double v_jp1 = in_fluid(i, jp1, k) ? tmp[idx(i,jp1,k)] : tmp[idx(i,j,k)];
                        field.x[i][j][k] = tmp[idx(i,j,k)]
                            + coeff * (v_jm1 - 2.0 * tmp[idx(i,j,k)] + v_jp1);
                    }
                }
            }
        }

        // --- along i (vertical, clamped at fluid floor and top) ---
        if (along_i) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    // clamp to [0, im-1]: guards against i_bathymetry < 0 when depth > L_hyd
                    const int i_surf = i_surface
                        ? std::max((*i_surface)[j][k], 0) : 0;
                    for (int i = i_surf; i < im; ++i) {
                        const int im1 = std::max(i - 1, i_surf);
                        const int ip1 = std::min(i + 1, im - 1);
                        field.x[i][j][k] = tmp[idx(i,j,k)]
                            + coeff * (tmp[idx(im1,j,k)]
                                     - 2.0 * tmp[idx(i,j,k)]
                                     + tmp[idx(ip1,j,k)]);
                    }
                }
            }
        }

    } // passes
}



// ============================================================================
// Latitude-dependent zonal (φ) "polar filter"
// ============================================================================
// Periodic 1-2-1 Shapiro step along k only, with a per-row pass count that grows
// toward the poles so the effective zonal resolution is held to that at the
// reference latitude.  See the header for the rationale (lat-lon polar CFL).
void AtomUtils::polar_zonal_filter(Array& field,
                                   const double* colat,
                                   const std::vector<std::vector<int>>* i_surface,
                                   double lat_filter_start_deg,
                                   int    max_passes,
                                   double pass_gain)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (max_passes <= 0 || !colat || pass_gain <= 0.0) return;

    // Reference metric factor sinθ_ref = cos(lat_filter_start) and a per-row pass
    // count passes[j] = clamp( round(sinθ_ref/sinθ_j − 1), 0, max_passes ).
    const double deg2rad = M_PI / 180.0;
    const double sin_ref = std::cos(lat_filter_start_deg * deg2rad);   // = sinθ at the reference latitude

    std::vector<int> passes(jm, 0);
    int global_max = 0;
    for (int j = 0; j < jm; ++j) {
        const double sin_j = std::max(std::sin(colat[j]), 1.0e-6);     // guard the pole rows (sinθ→0)
        if (sin_j >= sin_ref) continue;                                // equatorward of the reference latitude: untouched
        const double ratio = sin_ref / sin_j;                          // >1 poleward of the reference latitude
        int p = (int)std::lround(pass_gain * (ratio * ratio - 1.0));   // passes ∝ ratio² (√n smoothing length, see header), scaled by pass_gain
        if (p < 0) p = 0;
        if (p > max_passes) p = max_passes;
        passes[j] = p;
        if (p > global_max) global_max = p;
    }
    if (global_max == 0) return;

    auto in_fluid = [&](int i, int j, int k) -> bool {
        if (!i_surface) return true;
        return i >= std::max((*i_surface)[j][k], 0);
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };

    const double coeff = 0.25;   // full 1-2-1 Shapiro step (strength = 1)

    for (int pass = 0; pass < global_max; ++pass) {

        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < im; ++i) {
            for (int j = 0; j < jm; ++j) {
                if (passes[j] <= pass) continue;                        // this row has reached its pass quota
                for (int k = 0; k < km; ++k) {
                    if (!in_fluid(i, j, k)) continue;
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    // no-flux at solid neighbours: substitute current cell value
                    const double v_km1 = in_fluid(i, j, km1) ? tmp[idx(i,j,km1)] : tmp[idx(i,j,k)];
                    const double v_kp1 = in_fluid(i, j, kp1) ? tmp[idx(i,j,kp1)] : tmp[idx(i,j,k)];
                    field.x[i][j][k] = tmp[idx(i,j,k)]
                        + coeff * (v_km1 - 2.0 * tmp[idx(i,j,k)] + v_kp1);
                }
            }
        }
    } // passes
}



// ============================================================================
// Orographic Shapiro filter (steep-terrain CFL damping at any latitude)
// ============================================================================
// See the header for rationale. Flags columns whose surface index jumps by
// >= steep_threshold cells to a horizontal neighbour, then applies `passes`
// periodic 1-2-1 Shapiro sweeps along k and j to the air cells within
// n_layers_above of the surface there, with no-flux at solid neighbours.
void AtomUtils::orographic_shapiro_filter(Array& field,
                                          const std::vector<std::vector<int>>& i_surface,
                                          int    steep_threshold,
                                          int    n_layers_above,
                                          int    passes,
                                          double strength)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (passes <= 0 || strength <= 0.0) return;

    // Flag steep columns: |Δ surface-index| to any of the 4 horizontal neighbours
    // (k periodic, j clamped at the poles) >= steep_threshold cells.
    std::vector<unsigned char> steep(jm * km, 0);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; ++j) {
        for (int k = 0; k < km; ++k) {
            const int s   = std::max(i_surface[j][k], 0);
            const int jm1 = (j > 0)      ? j - 1 : 0;
            const int jp1 = (j < jm - 1) ? j + 1 : jm - 1;
            const int km1 = (k - 1 + km) % km;
            const int kp1 = (k + 1)      % km;
            int d = std::abs(s - std::max(i_surface[jm1][k], 0));
            d = std::max(d, std::abs(s - std::max(i_surface[jp1][k], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][km1], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][kp1], 0)));
            if (d >= steep_threshold) steep[j * km + k] = 1;
        }
    }

    auto in_fluid = [&](int i, int j, int k) {
        return i >= std::max(i_surface[j][k], 0);
    };
    // Cells we actually update: steep column, fluid, within n_layers_above of surface.
    auto target = [&](int i, int j, int k) {
        return steep[j * km + k] && in_fluid(i, j, k)
            && i < std::max(i_surface[j][k], 0) + n_layers_above;
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };
    const double coeff = 0.25 * strength;

    for (int pass = 0; pass < passes; ++pass) {

        // --- zonal (k, periodic) sweep ---
        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < im; ++i) {
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    if (!target(i, j, k)) continue;
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    const double c  = tmp[idx(i,j,k)];
                    const double vm = in_fluid(i, j, km1) ? tmp[idx(i,j,km1)] : c;
                    const double vp = in_fluid(i, j, kp1) ? tmp[idx(i,j,kp1)] : c;
                    field.x[i][j][k] = c + coeff * (vm - 2.0 * c + vp);
                }
            }
        }

        // --- meridional (j) sweep ---
        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < im; ++i) {
            for (int j = 1; j < jm - 1; ++j) {
                for (int k = 0; k < km; ++k) {
                    if (!target(i, j, k)) continue;
                    const double c  = tmp[idx(i,j,k)];
                    const double vm = in_fluid(i, j-1, k) ? tmp[idx(i,j-1,k)] : c;
                    const double vp = in_fluid(i, j+1, k) ? tmp[idx(i,j+1,k)] : c;
                    field.x[i][j][k] = c + coeff * (vm - 2.0 * c + vp);
                }
            }
        }

    } // passes
}



// ============================================================================
// Soft steep-massif field smoothing (velocity + pressure over steep contours)
// ============================================================================
// Companion to the init-time massif TERRAIN smoothing: even with the steep ramparts cut,
// the residual steepest columns keep a locally unphysical field aloft (user-observed over
// the Himalaya: saturated p_dyn block + bad u/v/w), while smooth-contour terrain like the
// South Pole behaves well. This applies a GENTLE horizontal (k then j) 1-2-1 Shapiro to a
// field ONLY over HIGH+STEEP contour columns (same mask as the terrain smoothing:
// i_surface ≥ high_thresh AND |Δi_surface to a neighbour| ≥ steep_thresh), reaching
// n_layers_deep cells up the column so it covers the troposphere above the steep mountain —
// NOT just the boundary layer. The mask is narrow (steep mountain interiors only), so plains,
// coasts, broad high terrain (South Pole) and the open ocean / i=0 vortices are untouched:
// it smooths the unphysical near-orography field without touching the expected real flow.
// Apply to u,v,w and p_dyn. No-flux (centre substitution) at land/sub-surface neighbours.
void AtomUtils::steep_massif_field_smoothing(Array& field,
                                             const std::vector<std::vector<int>>& i_surface,
                                             int    high_thresh,
                                             int    steep_thresh,
                                             int    n_layers_deep,
                                             int    passes,
                                             double strength)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (passes <= 0 || strength <= 0.0 || n_layers_deep <= 0 || jm < 3 || km < 3) return;

    // High+steep contour mask (identical to the init-time terrain massif smoothing).
    std::vector<unsigned char> mask(jm * km, 0);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; ++j) {
        for (int k = 0; k < km; ++k) {
            const int s = std::max(i_surface[j][k], 0);
            if (s < high_thresh) continue;
            const int jm1 = (j > 0)      ? j - 1 : 0;
            const int jp1 = (j < jm - 1) ? j + 1 : jm - 1;
            const int km1 = (k - 1 + km) % km;
            const int kp1 = (k + 1)      % km;
            int d = std::abs(s - std::max(i_surface[jm1][k], 0));
            d = std::max(d, std::abs(s - std::max(i_surface[jp1][k], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][km1], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][kp1], 0)));
            if (d >= steep_thresh) mask[j * km + k] = 1;
        }
    }

    auto in_fluid = [&](int i, int j, int k) {
        return i >= std::max(i_surface[j][k], 0);
    };
    auto target = [&](int i, int j, int k) {
        return mask[j * km + k] && in_fluid(i, j, k)
            && i < std::max(i_surface[j][k], 0) + n_layers_deep;
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };
    const double coeff = 0.25 * strength;

    for (int pass = 0; pass < passes; ++pass) {

        // --- zonal (k, periodic) sweep ---
        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < im; ++i) {
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    if (!target(i, j, k)) continue;
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    const double c  = tmp[idx(i,j,k)];
                    const double vm = in_fluid(i, j, km1) ? tmp[idx(i,j,km1)] : c;
                    const double vp = in_fluid(i, j, kp1) ? tmp[idx(i,j,kp1)] : c;
                    field.x[i][j][k] = c + coeff * (vm - 2.0 * c + vp);
                }
            }
        }

        // --- meridional (j) sweep ---
        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < im; ++i) {
            for (int j = 1; j < jm - 1; ++j) {
                for (int k = 0; k < km; ++k) {
                    if (!target(i, j, k)) continue;
                    const double c  = tmp[idx(i,j,k)];
                    const double vm = in_fluid(i, j-1, k) ? tmp[idx(i,j-1,k)] : c;
                    const double vp = in_fluid(i, j+1, k) ? tmp[idx(i,j+1,k)] : c;
                    field.x[i][j][k] = c + coeff * (vm - 2.0 * c + vp);
                }
            }
        }
    } // passes
}



// ============================================================================
// Gentle global radial (i) Shapiro de-checkerboarding
// ============================================================================
// The steep-orography velocity blow-up is a 2Δ grid-scale CHECKERBOARD whose purest
// component is on the RADIAL axis (the wmode diagnostic measured r=1.000 in i at the
// growing cell: the field is huge at one level and ~0/opposite-sign above and below).
// No other filter touches i — polar_zonal_filter is k only, orographic_shapiro_filter is
// j+k AND only at steep-flagged columns, which do NOT include the near-surface cell where
// the checkerboard actually grows (51°N/11°E, north of the steep Alps). So the radial
// component grows undamped → the iter-250 runaway. This applies `passes` 1-2-1 radial
// sweeps to ALL fluid cells (not gated by steepness): a single 1-2-1 (coeff 0.25·strength)
// annihilates a level-to-level 2Δ oscillation in one pass while barely touching smooth
// vertical shear (large scales are reduced only ~(Δ/λ)² per pass), so it is safe to apply
// globally every iteration. No-flux (current-cell substitution) at the solid surface
// below; i=0 and i=im-1 are left to the radial BCs.
void AtomUtils::radial_shapiro_filter(Array& field,
                                      const std::vector<std::vector<int>>& i_surface,
                                      int    passes,
                                      double strength)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (passes <= 0 || strength <= 0.0 || im < 3) return;

    auto in_fluid = [&](int i, int j, int k) {
        return i >= std::max(i_surface[j][k], 0);
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };
    const double coeff = 0.25 * strength;

    for (int pass = 0; pass < passes; ++pass) {

        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < im - 1; ++i) {
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    if (!in_fluid(i, j, k)) continue;
                    const double c  = tmp[idx(i,j,k)];
                    const double vm = in_fluid(i-1, j, k) ? tmp[idx(i-1,j,k)] : c;
                    const double vp = in_fluid(i+1, j, k) ? tmp[idx(i+1,j,k)] : c;
                    field.x[i][j][k] = c + coeff * (vm - 2.0 * c + vp);
                }
            }
        }
    } // passes
}



// ============================================================================
// Higher-order (4th-order Shapiro) radial de-checkerboarding
// ============================================================================
// Response 1 - sin^4(kΔ/2): removes the 2Δ grid mode completely while preserving
// resolved vertical shear far better than the 1-2-1 (the 1-2-1 run every iter
// homogenizes the jet's thermal-wind shear and the Hadley/Ferrel two-branch v
// profile, spinning the circulation down despite steady baroclinic forcing).
// Interior 5-point stencil  f - (1/16)(f_{i+2} -4 f_{i+1} +6 f_i -4 f_{i-1} + f_{i-2});
// where i±2 is unavailable (next to the surface or the lid) it falls back to the
// gentle 1-2-1. Solid cells below are substituted with the current value (no-flux).
void AtomUtils::radial_shapiro_filter_ho(Array& field,
                                         const std::vector<std::vector<int>>& i_surface,
                                         int    passes,
                                         double strength)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (passes <= 0 || strength <= 0.0 || im < 3) return;

    auto in_fluid = [&](int i, int j, int k) {
        return i >= std::max(i_surface[j][k], 0);
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };
    const double c4 = strength / 16.0;   // 4th-order Shapiro coefficient
    const double c2 = strength * 0.25;   // 1-2-1 fallback coefficient

    for (int pass = 0; pass < passes; ++pass) {

        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < im - 1; ++i) {
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    if (!in_fluid(i, j, k)) continue;
                    const double c = tmp[idx(i,j,k)];
                    if (!in_fluid(i-1, j, k) || !in_fluid(i+1, j, k)) continue;  // 1-cell-thin: leave it
                    const bool fm2 = (i - 2 >= 0)      && in_fluid(i-2, j, k);
                    const bool fp2 = (i + 2 <= im - 1) && in_fluid(i+2, j, k);
                    if (fm2 && fp2) {
                        const double d4 = tmp[idx(i+2,j,k)] - 4.0 * tmp[idx(i+1,j,k)]
                                        + 6.0 * c - 4.0 * tmp[idx(i-1,j,k)] + tmp[idx(i-2,j,k)];
                        field.x[i][j][k] = c - c4 * d4;          // 4th-order: kills 2Δ, keeps shear
                    } else {
                        const double vm = tmp[idx(i-1,j,k)];
                        const double vp = tmp[idx(i+1,j,k)];
                        field.x[i][j][k] = c + c2 * (vm - 2.0 * c + vp);  // boundary fallback
                    }
                }
            }
        }
    } // passes
}



// ============================================================================
// Orographic radial (i) scalar de-checkerboarding (extremum-gated)
// ============================================================================
// The NZ-Alps late blow-up is seeded by a single-cell 2Δz vertical TEMPERATURE
// checkerboard at the orographic condensation level (cold cell at i_surface+7 ≈ 937 m
// over the steep slope-capped Alps, 44°S/169°E), present from the dry phase
// (~iter150) long before moist physics; it later drives a moisture + vertical-velocity
// runaway via latent heat + buoyancy ([[project_nz_alps_t_checkerboard]]). minmod limits
// only ADVECTION, and the global radial Shapiro is applied to velocity only — nothing
// damps a radial 2Δz mode in the SCALAR fields. This applies a radial 1-2-1 to scalars
// but, to avoid the global-field perturbation that got the 2026-06-08 horizontal scalar
// orographic_shapiro_filter reverted (it tipped the 62°N seam), it acts ONLY at:
//   (1) steep/sloped columns (|Δ surface-index| ≥ steep_threshold), the orographic hot-spots,
//   (2) the first n_layers_above fluid cells, where the orographic mode lives, and
//   (3) cells that are a vertical local extremum ((vm−c)(vp−c) > 0) — the de-checkerboarding
//       analogue of minmod, so a monotonic lapse rate / inversion gets ZERO smoothing and
//       smooth columns are left bit-unchanged.
// No-flux (current-cell substitution) into solid cells below; i=0/i=im−1 left to the BCs.
void AtomUtils::orographic_radial_shapiro_filter(Array& field,
                                                 const std::vector<std::vector<int>>& i_surface,
                                                 int    steep_threshold,
                                                 int    n_layers_above,
                                                 int    passes,
                                                 double strength)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (passes <= 0 || strength <= 0.0 || im < 3) return;

    // Flag steep columns: |Δ surface-index| to any of the 4 horizontal neighbours >= threshold.
    std::vector<unsigned char> steep(jm * km, 0);
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < jm; ++j) {
        for (int k = 0; k < km; ++k) {
            const int s   = std::max(i_surface[j][k], 0);
            const int jm1 = (j > 0)      ? j - 1 : 0;
            const int jp1 = (j < jm - 1) ? j + 1 : jm - 1;
            const int km1 = (k - 1 + km) % km;
            const int kp1 = (k + 1)      % km;
            int d = std::abs(s - std::max(i_surface[jm1][k], 0));
            d = std::max(d, std::abs(s - std::max(i_surface[jp1][k], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][km1], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][kp1], 0)));
            if (d >= steep_threshold) steep[j * km + k] = 1;
        }
    }

    auto in_fluid = [&](int i, int j, int k) {
        return i >= std::max(i_surface[j][k], 0);
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };
    const double coeff = 0.25 * strength;

    for (int pass = 0; pass < passes; ++pass) {

        for (int i = 0; i < im; ++i)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[idx(i,j,k)] = field.x[i][j][k];

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < im - 1; ++i) {
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    if (!steep[j * km + k]) continue;
                    const int s = std::max(i_surface[j][k], 0);
                    if (i < s || i >= s + n_layers_above) continue;
                    const double c  = tmp[idx(i,j,k)];
                    const double vm = in_fluid(i-1, j, k) ? tmp[idx(i-1,j,k)] : c;
                    const double vp = in_fluid(i+1, j, k) ? tmp[idx(i+1,j,k)] : c;
                    if ((vm - c) * (vp - c) > 0.0)   // vertical local extremum → 2Δz checkerboard
                        field.x[i][j][k] = c + coeff * (vm - 2.0 * c + vp);
                }
            }
        }
    } // passes
}



// ============================================================================
// Near-surface coastal Rayleigh sponge
// ============================================================================
// Soft cap for the dry-seeded coastal velocity runaway (Gulf-of-Alaska mode).
// Flags any column whose i_surface neighbour-step >= 1 (every coast + slope),
// then at the lowest `layers` air cells, if |field| > threshold, multiplies by
// (1 - rate). Below threshold: untouched (preserves real coastal circulation).
void AtomUtils::coastal_velocity_sponge(Array& field,
                                        const std::vector<std::vector<int>>& i_surface,
                                        double threshold,
                                        double rate,
                                        int    layers)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (rate <= 0.0 || threshold < 0.0 || layers <= 0 || im < 2) return;

    std::vector<char> coastal(jm * km, 0);
    for (int j = 0; j < jm; ++j) {
        for (int k = 0; k < km; ++k) {
            const int s   = std::max(i_surface[j][k], 0);
            const int jm1 = std::max(j - 1, 0);
            const int jp1 = std::min(j + 1, jm - 1);
            const int km1 = (k == 0) ? km - 1 : k - 1;
            const int kp1 = (k == km - 1) ? 0 : k + 1;
            int d = std::abs(s - std::max(i_surface[jm1][k], 0));
            d = std::max(d, std::abs(s - std::max(i_surface[jp1][k], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][km1], 0)));
            d = std::max(d, std::abs(s - std::max(i_surface[j][kp1], 0)));
            if (d >= 1) coastal[j * km + k] = 1;
        }
    }

    // rate >= 1.0 → hard clamp at threshold (the soft Rayleigh form is multiplicative
    // and cannot bound the saturated coastal runaway where per-step RHS forcing >> cap;
    // a hard clamp acts as an unbreakable ceiling). rate < 1.0 → Rayleigh: v *= (1-rate)
    const bool hard_clamp = (rate >= 1.0);
    const double keep = hard_clamp ? 0.0 : (1.0 - rate);

    int n_coastal = 0;
    for (int j = 0; j < jm; ++j) for (int k = 0; k < km; ++k) if (coastal[j*km+k]) n_coastal++;
    long n_clamped = 0;
    double max_pre = 0.0;

    // Alaska diagnostic: 60°N/151°W ≈ (j=30, k=209) — but the print location varies
    // with sign conventions, so sample a 5×5 area and report the cell with max |v|.
    int j_ak_lo = 28, j_ak_hi = 32;
    int k_ak_lo = 207, k_ak_hi = 211;
    int alaska_flagged = 0;
    for (int j = j_ak_lo; j <= j_ak_hi; ++j)
        for (int k = k_ak_lo; k <= k_ak_hi; ++k)
            if (coastal[j * km + k]) alaska_flagged++;
    double alaska_max_pre = 0.0;
    int alaska_jmax = -1, alaska_kmax = -1, alaska_imax = -1, alaska_topo = -1;
    for (int j = j_ak_lo; j <= j_ak_hi; ++j) {
        for (int k = k_ak_lo; k <= k_ak_hi; ++k) {
            const int s = std::max(i_surface[j][k], 0);
            for (int i = 0; i < im; ++i) {
                const double v = std::abs(field.x[i][j][k]);
                if (v > alaska_max_pre) {
                    alaska_max_pre = v;
                    alaska_jmax = j; alaska_kmax = k; alaska_imax = i;
                    alaska_topo = s;
                }
            }
        }
    }

    #pragma omp parallel for collapse(2) schedule(static) reduction(+:n_clamped) reduction(max:max_pre)
    for (int j = 0; j < jm; ++j) {
        for (int k = 0; k < km; ++k) {
            if (!coastal[j * km + k]) continue;
            const int s = std::max(i_surface[j][k], 0);
            const int i_end = std::min(im, s + layers);
            for (int i = s; i < i_end; ++i) {
                const double v = field.x[i][j][k];
                const double av = std::abs(v);
                if (av > max_pre) max_pre = av;
                if (av > threshold) {
                    ++n_clamped;
                    if (hard_clamp) {
                        field.x[i][j][k] = (v > 0.0) ? threshold : -threshold;
                    } else {
                        field.x[i][j][k] = v * keep;
                    }
                }
            }
        }
    }
    // GLOBAL scan for the cell with max |v| at i=5 (the print's "298m" level).
    // Reports whether the global hot-cell is flagged by the coastal mask.
    double i5_max = 0.0;
    int    i5_j = -1, i5_k = -1, i5_topo = -1, i5_flag = -1;
    for (int j = 0; j < jm; ++j) {
        for (int k = 0; k < km; ++k) {
            const int s = std::max(i_surface[j][k], 0);
            if (s > 5) continue;
            const double v = std::abs(field.x[5][j][k]);
            if (v > i5_max) {
                i5_max = v;
                i5_j = j; i5_k = k; i5_topo = s;
                i5_flag = coastal[j * km + k];
            }
        }
    }

    std::cout << "[coastal_sponge] flagged=" << n_coastal
              << " clamped=" << n_clamped
              << " max_pre=" << max_pre
              << " AK_max=" << alaska_max_pre
              << " i5_max=" << i5_max
              << " @(j=" << i5_j << ",k=" << i5_k << ",topo=" << i5_topo << ",flag=" << i5_flag << ")"
              << std::endl;
}



// ============================================================================
// Upper-troposphere Rayleigh sponge
// ============================================================================
// Damps zonal anomalies (departures from the zonal mean of `field` at each (i,j))
// in the upper levels i ≥ i_start. Standard GCM upper-atmosphere remedy for the
// lat-lon CFL near the meridian seam at high altitude.
void AtomUtils::upper_rayleigh_sponge(Array& field, int i_start, double alpha_max)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    if (i_start >= im - 1 || alpha_max <= 0.0 || km <= 0) return;
    const double i_span = static_cast<double>(im - 1 - i_start);
    if (i_span <= 0.0) return;
    const double inv_km = 1.0 / static_cast<double>(km);

    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = i_start; i < im; i++) {
        for (int j = 0; j < jm; j++) {
            // Zonal mean at (i, j)
            double sum = 0.0;
            for (int k = 0; k < km; k++) sum += field.x[i][j][k];
            const double zonal_mean = sum * inv_km;

            // Linear ramp: alpha=0 at i_start, alpha=alpha_max at top
            const double alpha = alpha_max * (i - i_start) / i_span;
            const double keep  = 1.0 - alpha;

            for (int k = 0; k < km; k++) {
                field.x[i][j][k] = keep * field.x[i][j][k] + alpha * zonal_mean;
            }
        }
    }
}


// ============================================================================
// Extreme-peak removal filter
// ============================================================================
// Along each enabled axis the two direct neighbours define a local background
// (mean) and variability scale (half-spread = |f[n+1]-f[n-1]|/2).
// Cells that deviate from the background by more than threshold*half_spread
// are replaced by the background.  Smooth gradients are untouched.
void AtomUtils::remove_peaks(Array& field,
                              const std::vector<std::vector<int>>* i_surface,
                              bool along_i, bool along_j, bool along_k,
                              double threshold,
                              int passes)
{
    const int im = field.im;
    const int jm = field.jm;
    const int km = field.km;

    auto in_fluid = [&](int i, int j, int k) -> bool {
        if (!i_surface) return true;
        return i >= std::max((*i_surface)[j][k], 0);
    };

    std::vector<double> tmp(im * jm * km);
    auto idx = [&](int i, int j, int k) { return i * jm * km + j * km + k; };

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        const int km1 = (k - 1 + km) % km;
                        const int kp1 = (k + 1)      % km;
                        const double v   = tmp[idx(i,j,k)];
                        const double v_m = in_fluid(i, j, km1) ? tmp[idx(i,j,km1)] : v;
                        const double v_p = in_fluid(i, j, kp1) ? tmp[idx(i,j,kp1)] : v;
                        const double mean        = 0.5 * (v_m + v_p);
                        const double half_spread = 0.5 * std::abs(v_p - v_m);
                        if (std::abs(v - mean) > threshold * half_spread)
                            field.x[i][j][k] = mean;
                    }
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int i = 0; i < im; ++i) {
                for (int j = 0; j < jm; ++j) {
                    const int jm1 = std::max(j - 1, 0);
                    const int jp1 = std::min(j + 1, jm - 1);
                    for (int k = 0; k < km; ++k) {
                        if (!in_fluid(i, j, k)) continue;
                        const double v   = tmp[idx(i,j,k)];
                        const double v_m = in_fluid(i, jm1, k) ? tmp[idx(i,jm1,k)] : v;
                        const double v_p = in_fluid(i, jp1, k) ? tmp[idx(i,jp1,k)] : v;
                        const double mean        = 0.5 * (v_m + v_p);
                        const double half_spread = 0.5 * std::abs(v_p - v_m);
                        if (std::abs(v - mean) > threshold * half_spread)
                            field.x[i][j][k] = mean;
                    }
                }
            }
        }

        // --- along i (vertical, clamped at fluid floor and top) ---
        if (along_i) {
            for (int i = 0; i < im; ++i)
                for (int j = 0; j < jm; ++j)
                    for (int k = 0; k < km; ++k)
                        tmp[idx(i,j,k)] = field.x[i][j][k];

            #pragma omp parallel for collapse(2) schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    const int i_surf = i_surface
                        ? std::max((*i_surface)[j][k], 0) : 0;
                    for (int i = i_surf; i < im; ++i) {
                        const int im1 = std::max(i - 1, i_surf);
                        const int ip1 = std::min(i + 1, im - 1);
                        const double v   = tmp[idx(i,j,k)];
                        const double v_m = tmp[idx(im1,j,k)];
                        const double v_p = tmp[idx(ip1,j,k)];
                        const double mean        = 0.5 * (v_m + v_p);
                        const double half_spread = 0.5 * std::abs(v_p - v_m);
                        if (std::abs(v - mean) > threshold * half_spread)
                            field.x[i][j][k] = mean;
                    }
                }
            }
        }

    } // passes
}

// ============================================================================
// Shapiro (1-2-1) wiggle-damping filter — Array_2D overload
// ============================================================================
void AtomUtils::damp_wiggles(Array_2D& field,
                              bool along_j, bool along_k,
                              double strength,
                              int passes)
{
    const int jm = field.jm;
    const int km = field.km;
    const double coeff = strength * 0.25;

    // Flat scratch buffer: tmp[j*km + k]
    std::vector<double> tmp(jm * km);

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    field.y[j][k] = tmp[j * km + k]
                        + coeff * (tmp[j * km + km1]
                                 - 2.0 * tmp[j * km + k]
                                 + tmp[j * km + kp1]);
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                const int jm1 = std::max(j - 1, 0);
                const int jp1 = std::min(j + 1, jm - 1);
                for (int k = 0; k < km; ++k) {
                    field.y[j][k] = tmp[j * km + k]
                        + coeff * (tmp[jm1 * km + k]
                                 - 2.0 * tmp[j   * km + k]
                                 + tmp[jp1 * km + k]);
                }
            }
        }

    } // passes
}

// ============================================================================
// Extreme-peak removal filter — Array_2D overload
// ============================================================================
void AtomUtils::remove_peaks(Array_2D& field,
                              bool along_j, bool along_k,
                              double threshold,
                              int passes)
{
    const int jm = field.jm;
    const int km = field.km;

    std::vector<double> tmp(jm * km);

    for (int pass = 0; pass < passes; ++pass) {

        // --- along k (longitude, periodic) ---
        if (along_k) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                for (int k = 0; k < km; ++k) {
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    const double v   = tmp[j * km + k];
                    const double v_m = tmp[j * km + km1];
                    const double v_p = tmp[j * km + kp1];
                    const double mean        = 0.5 * (v_m + v_p);
                    const double half_spread = 0.5 * std::abs(v_p - v_m);
                    if (std::abs(v - mean) > threshold * half_spread)
                        field.y[j][k] = mean;
                }
            }
        }

        // --- along j (latitude, clamped at poles) ---
        if (along_j) {
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j)
                for (int k = 0; k < km; ++k)
                    tmp[j * km + k] = field.y[j][k];

            #pragma omp parallel for schedule(static)
            for (int j = 0; j < jm; ++j) {
                const int jm1 = std::max(j - 1, 0);
                const int jp1 = std::min(j + 1, jm - 1);
                for (int k = 0; k < km; ++k) {
                    const double v   = tmp[j   * km + k];
                    const double v_m = tmp[jm1 * km + k];
                    const double v_p = tmp[jp1 * km + k];
                    const double mean        = 0.5 * (v_m + v_p);
                    const double half_spread = 0.5 * std::abs(v_p - v_m);
                    if (std::abs(v - mean) > threshold * half_spread)
                        field.y[j][k] = mean;
                }
            }
        }

    } // passes
}

