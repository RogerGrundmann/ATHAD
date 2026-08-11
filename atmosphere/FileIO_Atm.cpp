#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <dirent.h>
#include <sstream>

#include "cAtmosphereModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;
 

void cAtmosphereModel::read_Atmosphere_Surface_Data(int Ma){

cout << endl << endl << endl << "      AGCM: read_Atmosphere_Surface_Data ......................." << endl;

    auto begin = std::chrono::high_resolution_clock::now();


    if(!has_printed_welcome_msg)  print_welcome_msg();


// loading topography
    bathymetry_name = std::to_string(Ma) + BathymetrySuffix;

    std::cout << endl << "      topography given by the x-y-z data set:    " 
        << bathymetry_name.c_str() << endl;

    init_topography(bathymetry_path + "/" + bathymetry_name);
    { auto dot = bathymetry_name.rfind('.'); if(dot != string::npos) bathymetry_name.erase(dot); }




// temperature file loading for around 500 Mio years
    #pragma omp parallel sections
    {
        #pragma omp section
        { if(!is_global_temperature_curve_loaded()) load_global_temperature_curve(); }// loading of global temperature files from Scotese et al., ESR 2021 in °C
        #pragma omp section
        { if(!is_equat_temperature_curve_loaded()) load_equat_temperature_curve(); }  // loading of equatorial temperature files from Scotese et al., ESR 2021 in °C
        #pragma omp section
        { if(!is_pole_temperature_curve_loaded()) load_pole_temperature_curve(); }    // loading of polar temperature files from Scotese et al., ESR 2021 in °C
    }


// preparing temperature, precipitation and velocity data file constructed by reconstruct_atom_data.py
    string Name_SurfaceTemperature_File  = temperature_file;            // if  Ma > 0, temperature file name
    string Name_SurfaceNASATemperature_File  = temperature_file;        // if  Ma = 0, NASA temperature file name
    string Name_SurfaceNASAPrecipitation_File = precipitation_file;     // precipitation file name

    if((*get_current_time() > 0)&&(use_earthbyte_reconstruction)){      // reading temperature, precipitation, v and w in case Ma > 0
        Name_SurfaceTemperature_File = output_path
            + std::to_string(Ma) + "Ma_Reconstructed_Temperature.xyz";  // reconstructed temperature in °C

        Name_SurfaceNASAPrecipitation_File = output_path
            + std::to_string(Ma) + "Ma_Reconstructed_Precipitation.xyz";// reconstructed precipitation in mm/d

        velocity_v_file = output_path + std::to_string(Ma)
            + "Ma_Reconstructed_wind_v.xyz";                            // reconstructed v-velocity in m/s

        velocity_w_file = output_path + std::to_string(Ma)
            + "Ma_Reconstructed_wind_w.xyz";                            // reconstructed w-velocity in m/s


        struct stat info;                                               // initiates the search of corrupted reconstruction files

        if(stat(output_path.c_str(), &info) != 0){                      // if output path does not exist mkdir is applied
             mkdir(output_path.c_str(), 0777);
        }


        if(stat(Name_SurfaceTemperature_File.c_str(), &info) != 0 ||    // if any file is corrupt an error message is posted,  0 == executed successfully
           stat(Name_SurfaceNASAPrecipitation_File.c_str(), &info) != 0 ||
           stat(velocity_v_file.c_str(), &info) != 0 ||
           stat(velocity_w_file.c_str(), &info) != 0){

            std::string cmd_str = "python " + reconstruction_script_path// python ../reconstruction/reconstruct_atom_data.py, runs the python code
                + " " + std::to_string(Ma - time_step)                  // sys.argv[1], preceding timestep
                + " " + std::to_string(Ma)                              // sys.argv[2], recent timestep
                + " " + output_path                                     // sys.argv[3], output path where to store reconstructed files
                + " " + BathymetrySuffix                                // sys.argv[4], BathymetrySuffix here 'Ma_smooth.xyz'
                + " atm";                                               // sys.argv[5], here Atmosphere code (atm or hyd) applied

            int ret = system(cmd_str.c_str());                          // executes an external operating system command stored in std::string cmd_str, here reconstruct_atom_data.py

            std::cout << " reconstruction script returned: "            // return value ret = 0 means successful ended, ret = -1 means the command processor for python cannot be created
                << ret << std::endl;
        } 
    }


    #pragma omp parallel sections
    {
        #pragma omp section
        {
            if((*get_current_time() > 0) && (use_earthbyte_reconstruction)){
                read_IC(Name_SurfaceTemperature_File, t.x[0], jm, km);
            }
        }
        #pragma omp section
        {
            read_IC(Name_SurfaceNASATemperature_File, temperature_NASA.y, jm, km);
            read_IC(Name_SurfaceNASAPrecipitation_File, precipitation_NASA.y, jm, km);

//            read_IC(velocity_v_file, velocity_v_NASA.y, jm, km);                // reconstructed v-velocity in m/s

//            read_IC(velocity_w_file, velocity_w_NASA.y, jm, km);                // reconstructed v-velocity in m/s
        }
    }


    cout << "      AGCM: read_Atmosphere_Surface_Data ended ................." << endl << endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for read_Atmosphere_Surface_Data\n", elapsed.count() * 1e-9);

}
/*
*
*/
void cAtmosphereModel::AtmosphereDataTransfer(const string &Name_Bathymetry_File){
    cout << endl << "      AGCM: AtmosphereDataTransfer   iter_n = " << iter_n << endl;

    string stem_t = Name_Bathymetry_File;
    { auto dot = stem_t.rfind('.'); if(dot != string::npos) stem_t = stem_t.substr(0, dot); }
    string base_t = output_path;
    if(!base_t.empty() && base_t.back() == '/') base_t.pop_back();
    string Name_Transfer_File = base_t + "/" + stem_t + "_Transfer_Atm.vwtp";

    // Build the buffer BEFORE touching the file, while checking every value that
    // would be written for NaN/Inf (is_finite_safe, since -ffast-math defeats
    // std::isfinite). If any is non-finite the file is NOT overwritten further
    // below, so the last clean version stays on disk for the hydrosphere to use.
    std::vector<std::string> line_buffer(jm * km);
    bool has_nonfinite = false;

    #pragma omp parallel for collapse(2) schedule(static) reduction(||: has_nonfinite)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            std::stringstream ss;
            ss.precision(4);
            ss.setf(ios::fixed);

            if(is_land(h, 0, j, k)){
                ss << "0.0000 0.0000 0.0000 0.0000 0.0000 0.0000";
            }else{
                double vv = v.x[0][j][k];                               //non-dimensional
                double wv = w.x[0][j][k];                               //non-dimensional
                double tv = t.x[0][j][k];                               //non-dimensional
                double pv = p_dyn.x[0][j][k];                           //non-dimensional
                double ev = Evaporation.y[j][k];                        //dimensional in [mm/d]
                double pr = Precipitation.x[0][j][k] * 86400.0;         //dimensional in [mm/d]

                if(!is_finite_safe(vv) || !is_finite_safe(wv) || !is_finite_safe(tv)
                || !is_finite_safe(pv) || !is_finite_safe(ev) || !is_finite_safe(pr))
                    has_nonfinite = true;

                ss << vv << " " << wv << " " << tv << " "
                   << pv << " " << ev << " " << pr;
            }

            // write to buffer (index must be unique per thread)
            line_buffer[j * km + k] = ss.str();
        }
    }

    // Guard: a NaN/Inf in the surface fields means this transfer would corrupt the
    // file the hydrosphere reads. Skip the write and keep the previous clean file.
    if(has_nonfinite){
        cout << "      AGCM: AtmosphereDataTransfer SKIPPED at iter_n = " << iter_n
             << " — non-finite (NaN/Inf) surface values present; previous transfer file left intact" << endl;
        return;
    }

    // Write the finished buffer to an iter-stamped snapshot
    //   <stem>_Transfer_Atm_<iter_n>.vwtp  — retained (no clobber), so the hydrosphere can be
    //   driven by a chosen atmosphere iteration, or by the LATEST one wherever the atmosphere
    //   run happened to stop (see HydrosphereDataTransfer / atm_transfer_iter, which auto-selects
    //   the highest-iter snapshot by default). Mirrors the atm_restart_<iter>.bin cadence/naming.
    //   Kept as TEXT (surface-only, ~few MB, human-inspectable BC file). No fixed-name "latest"
    //   copy is written any more — the hydrosphere discovers the newest snapshot itself.
    const string stamped_name = Name_Transfer_File.substr(0, Name_Transfer_File.size() - 5)  // strip ".vwtp"
                              + "_" + std::to_string(iter_n) + ".vwtp";
    {
        ofstream Transfer_File(stamped_name);
        if(!Transfer_File.is_open()){
            cerr << "ERROR: could not open transfer file: " << stamped_name << "\n";
            abort();
        }
        // headline carrying the iteration this transfer corresponds to; the hydrosphere
        // reader (HydrosphereDataTransfer) consumes this single line before the data rows.
        Transfer_File << "# iter_n = " << iter_n << "\n";
        for(int i = 0; i < jm * km; i++){
            Transfer_File << line_buffer[i] << "\n";
        }
        Transfer_File.close();
    }
    cout << "      AGCM: AtmosphereDataTransfer ended (wrote " << stamped_name << ")" << endl;
}
/*
*
*/
// Reverse coupling channel (hydrosphere -> atmosphere), the atmosphere half of the
// atm<->hyd Picard loop. Reads the ocean SST a previous hydrosphere run wrote
// (HydrosphereSSTTransfer -> <stem>_Transfer_Hyd_SST_<iter>.vwtp) and blends it into the
// atmospheric ocean surface temperature t.x[0], so the NEXT atmosphere solve sees an
// ocean-dynamics-consistent SST instead of the frozen Scotese parabola. Called at init,
// AFTER initTemperatureData sets t.x[0] and BEFORE run_3D_loop snapshots t_eq — so the
// blended SST also becomes the Held-Suarez relaxation target for free (see run_3D_loop).
//
// Safety (this must never wreck a normal one-way paleo run):
//   - OPT-IN: sst_coupling_alpha <= 0 returns immediately (no read) -> bit-identical to
//     the one-way chain, and to round 0 of a Picard loop (no SST file exists yet anyway).
//   - ocean cells only (land surface T is the atmosphere's own business);
//   - every value finite-checked; a non-finite entry leaves that cell untouched;
//   - SST hard-clamped to [-1.8, 40] C (non-dimensional) so a bad ocean cell cannot
//     inject an absurd surface temperature;
//   - OCEAN-MEAN ANCHORED: the ocean-mean of t.x[0] is restored after blending, so the
//     coupling moves SST STRUCTURE around without shifting the global ocean energy
//     (the Scotese ocean-mean, the one thing we trust, is preserved exactly).
// Under-relaxation across rounds is the caller's job (sst_coupling_alpha ~0.3-0.5).
void cAtmosphereModel::read_Hydrosphere_SST(int Ma){
    if(sst_coupling_alpha <= 0.0) return;                              // OFF: one-way chain / Picard round 0

    cout << endl << "      AGCM: read_Hydrosphere_SST   alpha = " << sst_coupling_alpha << endl;

    string stem_t = bathymetry_name;                                  // already extension-stripped at init
    { auto dot = stem_t.rfind('.'); if(dot != string::npos) stem_t = stem_t.substr(0, dot); }
    string base_t = output_path;
    if(!base_t.empty() && base_t.back() == '/') base_t.pop_back();

    // Select which hydrosphere SST snapshot to read (mirrors HydrosphereDataTransfer's
    // atm_transfer_iter logic): hyd_sst_iter >= 0 reads exactly that iteration; < 0 reads
    // the highest-iter snapshot present (i.e. wherever the hydrosphere run stopped).
    const string prefix = stem_t + "_Transfer_Hyd_SST_";
    string Name_SST_File;
    if(hyd_sst_iter >= 0){
        Name_SST_File = base_t + "/" + prefix + std::to_string(hyd_sst_iter) + ".vwtp";
    }else{
        int best = -1;
        if(DIR *dir = opendir(base_t.c_str())){
            while(struct dirent *ent = readdir(dir)){
                const string fn = ent->d_name;
                if(fn.size() > prefix.size() + 5
                && fn.compare(0, prefix.size(), prefix) == 0
                && fn.compare(fn.size() - 5, 5, ".vwtp") == 0){
                    const string mid = fn.substr(prefix.size(), fn.size() - prefix.size() - 5);
                    if(!mid.empty() && mid.find_first_not_of("0123456789") == string::npos)
                        best = std::max(best, std::stoi(mid));
                }
            }
            closedir(dir);
        }
        if(best < 0){
            cout << "      AGCM: no hydrosphere SST snapshot (" << prefix
                 << "<iter>.vwtp) in " << base_t << " — skipping SST coupling (round 0?)" << endl;
            return;
        }
        Name_SST_File = base_t + "/" + prefix + std::to_string(best) + ".vwtp";
    }

    ifstream SST_File(Name_SST_File);
    if(!SST_File.is_open()){
        cout << "      AGCM: hydrosphere SST file not found, skipping SST coupling: "
             << Name_SST_File << endl;
        return;
    }
    cout << "      AGCM: hydrosphere SST file = " << Name_SST_File
         << (hyd_sst_iter >= 0 ? "  (hyd_sst_iter)" : "  (latest snapshot)") << endl;

    string line;
    if(getline(SST_File, line)){                                      // optional "# iter_n = ..." headline
        if(!line.empty() && line[0] == '#')
            cout << "      AGCM: SST file headline: " << line << endl;
        else { SST_File.clear(); SST_File.seekg(0, ios::beg); }
    }

    // Non-dimensional SST clamp band: T_nd = T_K / t_0, so [-1.8, 40] C maps as below.
    const double sst_lo = (t_0 - 1.8) / t_0;
    const double sst_hi = (t_0 + 40.0) / t_0;

    // Pass 1: ocean-mean of t.x[0] BEFORE blending (the invariant to preserve).
    double sum_before = 0.0; long n_ocean = 0;
    for(int j = 0; j < jm; j++)
        for(int k = 0; k < km; k++)
            if(!is_land(h, 0, j, k)){ sum_before += t.x[0][j][k]; ++n_ocean; }

    // Pass 2: read + blend over ocean cells.
    bool truncated = false;
    for(int j = 0; j < jm && !truncated; j++){
        for(int k = 0; k < km; k++){
            if(!getline(SST_File, line)){
                cerr << "WARNING: read_Hydrosphere_SST: file truncated at j="
                     << j << " k=" << k << "\n";
                truncated = true; break;
            }
            if(is_land(h, 0, j, k)) continue;                         // ocean only
            double sst = 0.0; { istringstream ss(line); ss >> sst; }
            if(!is_finite_safe(sst)) continue;                        // leave cell untouched
            sst = std::min(std::max(sst, sst_lo), sst_hi);            // clamp to [-1.8,40] C
            t.x[0][j][k] = (1.0 - sst_coupling_alpha) * t.x[0][j][k]
                         + sst_coupling_alpha * sst;
        }
    }
    SST_File.close();

    // Pass 3: ocean-mean anchor — shift every ocean cell so the ocean-mean is unchanged.
    if(n_ocean > 0){
        double sum_after = 0.0;
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                if(!is_land(h, 0, j, k)) sum_after += t.x[0][j][k];
        const double corr = (sum_before - sum_after) / static_cast<double>(n_ocean);
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++)
                if(!is_land(h, 0, j, k)) t.x[0][j][k] += corr;
        cout << "      AGCM: SST coupling applied over " << n_ocean
             << " ocean cells (mean-anchor shift " << corr << " nd)" << endl;
    }
    cout << "      AGCM: read_Hydrosphere_SST ended" << endl;
}
/*
*
*/
void cAtmosphereModel::AtmospherePlotData(const string &Name_Bathymetry_File){
    cout << endl << "      AGCM: AtmospherePlotData   iter_n = " << iter_n << endl;

    string stem = Name_Bathymetry_File;
    { auto dot = stem.rfind('.'); if(dot != string::npos) stem = stem.substr(0, dot); }
    string base = output_path;
    if(!base.empty() && base.back() == '/') base.pop_back();
    string Name_PlotData_File = base + "/" + stem + "_PlotData_Atm.xyz";
    ofstream PlotData_File(Name_PlotData_File);

    if(!PlotData_File.is_open()){
        cerr << "ERROR: could not open PlotData file " << __FILE__ << " at line " << __LINE__ << "\n";
        abort();
    }

    PlotData_File.precision(4);
    PlotData_File.setf(ios::fixed);

    // Header schreiben (single line, kept on one line so parsers using
    // skiprows=1 / skip_header=1 stay correct). The actual iteration number
    // iter_n is appended for later inspection of which iteration produced this file.
    PlotData_File << "lons(deg), lats(deg), topography(m), v-velocity(m/s), w-velocity(m/s), "
                  << "velocity-mag(m/s), temperature(Celsius), water_vapour(g/kg), "
                  << "precipitation(mm/d), precipitable water(mm), pressure-static (hPa), "
                  << "temp_landscape (Celsius), p_stat_landscape (hPa)"
                  << "   ## iter_n = " << iter_n << endl;

    // Puffer für alle Zeilen (km * jm)
    std::vector<std::string> buffer(km * jm);

    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            // Berechnungen parallel durchführen
            double vel_v = v.x[0][j][k] * u_0;
            double vel_w = w.x[0][j][k] * u_0;
            double vel_mag = sqrt(vel_v * vel_v + vel_w * vel_w); // schneller als pow()

            std::stringstream ss;
            ss.precision(4);
            ss.setf(ios::fixed);

            ss << k << " " << 90-j 
               << " " << h.x[0][j][k] 
               << " " << vel_v
               << " " << vel_w
               << " " << vel_mag
               << " " << t.x[0][j][k] * t_0 - t_0
               << " " << c.x[0][j][k] * 1000.0
               << " " << Precipitation.x[0][j][k] * 86400.0
               << " " << precipitable_water.y[j][k] 
               << " " << p_stat.x[0][j][k]
               << " " << temp_landscape.y[j][k]
               << " " << p_stat_landscape.y[j][k];

            buffer[k * jm + j] = ss.str();
        }
    }

    // Serielles Schreiben des Puffers
    for(const auto& line : buffer) {
        PlotData_File << line << "\n";
    }

    PlotData_File.close();
    cout << "      AGCM: AtmospherePlotData ended" << endl;
}
/*
*
*/
void cAtmosphereModel::init_topography(const string &topo_filename){

    cout << endl << endl << endl << "      AGCM: init_topography" << endl;

    ifstream ifile(topo_filename);

    if(! ifile.is_open()){
        std::cerr << "ERROR: could not open Name_Bathymetry_File file: " 
            <<  topo_filename << std::endl;
        abort();
    }

    double lon, lat, height;

    for(int j = 0; j < jm && !ifile.eof(); j++){
        for(int k = 0; k < km && !ifile.eof(); k++){

            height = -999.0;                                            // in case the height is NaN

            ifile >> lon >> lat >> height;

            if(!(height > 0.0)){
                h.x[0][j][k] = Topography.y[j][k] = 0.0;
           }else{
                Topography.y[j][k] = height;

                for(int i = 0; i < im; i++){

                    if(height > get_layer_height(i)){
                        h.x[i][j][k] = 1.0;
                    }else{
                        i_topography[j][k] = i-1;
                        break;
                    }
                }
            }

            if(ifile.fail()){
                ifile.clear();
                std::string tmp;
                std::getline(ifile, tmp);
                logger() << "bad data in topography at: " << lon << " " 
                    << lat << " " << tmp << std::endl;
            }   

        }
    }

//  reduction and smoothing of peaks and needles in the topography
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < im; i++){
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){

                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j-1, k))
                     &&(is_air(h, i, j+1, k)))){
                    h.x[i][j][k] = 0.0;
                }
                if((is_land(h, i, j, k))
                     &&((is_air(h, i, j, k-1))
                     &&(is_air(h, i, j, k+1)))){
                    h.x[i][j][k] = 0.0;
                }

            }
        }
    }

    // Reset i_topography (and derived height arrays) to the OCEAN default BEFORE the
    // recompute below, so the index is authoritatively re-derived from the (now peak-
    // smoothed) h mask. Without this, a column whose land was FULLY cleared by the peak/
    // needle smoothing above keeps the stale surface index it got at read time — there is
    // no land→air transition left for the recompute to overwrite it — producing a "phantom
    // mountain" (i_topography>=1 over open water). Those phantom columns are treated as
    // fluid by the pressure solver (which keys on h) yet as terrain by the RHS/BCs (which
    // key on i_topography), injecting an anomalous coastal Poisson pressure that drove the
    // westerly-coast v/w velocity overshoot and the iter-357 NaN
    // ([[project-coastal-pgrad-pathology]]). A phase-by-phase scan measured 2881 such
    // columns appearing at the peak-smooth step and 353 surviving the recompute without
    // this reset; with it, 0 remain.
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            i_topography[j][k] = 0;
            i_landscape[j][k]  = get_layer_height(0);
            Landscape.y[j][k]  = get_layer_height(0);
        }
    }

    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 1; i < im; i++){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                if((is_air(h, i, j, k))&&(is_land(h, i-1, j, k))){
                    i_topography[j][k] = i;
                    i_landscape[j][k] = get_layer_height(i);
                    Landscape.y[j][k] = get_layer_height(i);
                }
            }
        }
    }

//  Orographic slope cap: limit the surface-index step between neighbouring LAND
//  columns to at most max_surface_jump cells. A near-vertical topographic wall
//  (the Himalaya/Tibet rampart jumps ~9 cells over one grid step) is the seed of the
//  cliff-driven pressure-divergence runaway — the Poisson solver spreads that step's
//  divergence source up the column and it blows up aloft (first NaN seen at ~12 km over
//  Tibet). Capping the slope removes the seed while PRESERVING peak height (we only
//  RAISE the low side, never cut peaks) and never converting ocean to land (ocean
//  columns have i_topography==0 and are skipped, so coastlines are untouched). The cliff
//  gets a broader, grid-resolvable apron so the flow deflects around a real mountain
//  instead of an unresolved wall. Pure-raise relaxation converges to a slope-<=cap cone
//  around each peak that stops where it meets terrain already within the cap.
    {
        const int max_surface_jump = 4;                      // max land-land surface-index step (tuning knob). NOTE (2026-06-09): tightening 4→2 was TESTED and made things WORSE — from-scratch it crashed at iter 338 (vs 483) and relocated the NaN to the Rockies (34°N/115°W). Global apron-broadening raises more land and seeds NEW steep edges elsewhere — too blunt a lever for the orographic-convergence seed. Reverted to 4. See [[project_upper_velocity_secular_growth]].
        std::vector<std::vector<int>> surf = i_topography;   // working copy: first-air index (0 over ocean)
        bool changed = true;
        for(int guard = 0; changed && guard < im; ++guard){
            changed = false;
            std::vector<std::vector<int>> prev = surf;
            for(int j = 0; j < jm; ++j){
                const int jm1 = (j > 0)      ? j - 1 : 0;
                const int jp1 = (j < jm - 1) ? j + 1 : jm - 1;
                for(int k = 0; k < km; ++k){
                    if(prev[j][k] == 0) continue;            // ocean column: never raise (keep coastline)
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    int ns = prev[j][k];                     // ocean neighbours (==0) give 0-cap < ns, so no-op
                    ns = std::max(ns, prev[jm1][k] - max_surface_jump);
                    ns = std::max(ns, prev[jp1][k] - max_surface_jump);
                    ns = std::max(ns, prev[j][km1] - max_surface_jump);
                    ns = std::max(ns, prev[j][kp1] - max_surface_jump);
                    if(ns > prev[j][k]){ surf[j][k] = ns; changed = true; }
                }
            }
        }
        // Rebuild the land mask and derived surface arrays for every raised column.
        for(int j = 0; j < jm; ++j){
            for(int k = 0; k < km; ++k){
                const int ns = surf[j][k];
                if(ns > i_topography[j][k]){
                    for(int i = i_topography[j][k]; i < ns; ++i) h.x[i][j][k] = 1.0;
                    i_topography[j][k] = ns;
                    i_landscape[j][k]  = get_layer_height(ns);
                    Landscape.y[j][k]  = get_layer_height(ns);
                    if(Topography.y[j][k] < get_layer_height(ns))
                        Topography.y[j][k] = get_layer_height(ns);
                }
            }
        }
    }

//  Targeted massif smoothing (2026-06-10, user-chosen root fix). The pure-raise slope cap
//  above PRESERVES peak height, so the steep HIGH massifs (Himalaya/Tibet/Tian-Shan/Pamir)
//  keep a near-vertical rampart even at slope==max_surface_jump. The base horizontal
//  convergence at that rampart is the divergence source that saturates p_dyn at the ±ceiling
//  (user saw an "equally high" p_dyn block + unnatural upwelling + awful v,w over the
//  Himalaya, while the broad-but-not-steep South Pole stays clean) and drives the upper-
//  troposphere velocity runaway / iter-493 NaN at the Pamir column
//  ([[project-upper-velocity-secular-growth]]). Here we additionally SMOOTH i_topography
//  (gentle 5-point 1-2-1, a few passes) ONLY over columns that are BOTH high AND steep —
//  this LOWERS the peak and broadens the apron, cutting the slope (hence the convergence)
//  below what the pure-raise cap alone can. Confined to the steep-high massif so plains,
//  coasts, broad high terrain and the ocean are untouched. h.x is rebuilt consistently for
//  BOTH raised and lowered columns to avoid the i_topography/h desync phantom-mountain bug
//  ([[project-coastal-pgrad-pathology]]).
    {
        const int massif_high_thresh  = 6;   // only high orography (Pamir≈7; Tibet/Himalaya higher)
        const int massif_steep_thresh = 3;   // only steep columns (|Δ surface-index| to a neighbour)
        const int smooth_passes       = 2;   // gentle

        std::vector<std::vector<char>> massif(jm, std::vector<char>(km, 0));
        for(int j = 0; j < jm; ++j){
            const int jm1 = (j > 0)      ? j - 1 : 0;
            const int jp1 = (j < jm - 1) ? j + 1 : jm - 1;
            for(int k = 0; k < km; ++k){
                const int s = i_topography[j][k];
                if(s < massif_high_thresh) continue;
                const int km1 = (k - 1 + km) % km;
                const int kp1 = (k + 1)      % km;
                int d = std::abs(s - i_topography[jm1][k]);
                d = std::max(d, std::abs(s - i_topography[jp1][k]));
                d = std::max(d, std::abs(s - i_topography[j][km1]));
                d = std::max(d, std::abs(s - i_topography[j][kp1]));
                if(d >= massif_steep_thresh) massif[j][k] = 1;
            }
        }

        std::vector<std::vector<int>> surf = i_topography;
        for(int pass = 0; pass < smooth_passes; ++pass){
            std::vector<std::vector<int>> prev = surf;
            for(int j = 0; j < jm; ++j){
                const int jm1 = (j > 0)      ? j - 1 : 0;
                const int jp1 = (j < jm - 1) ? j + 1 : jm - 1;
                for(int k = 0; k < km; ++k){
                    if(!massif[j][k]) continue;
                    const int km1 = (k - 1 + km) % km;
                    const int kp1 = (k + 1)      % km;
                    const int sum = 4 * prev[j][k]
                                  + prev[jm1][k] + prev[jp1][k]
                                  + prev[j][km1] + prev[j][kp1];
                    surf[j][k] = (sum + 4) / 8;   // rounded 1-2-1 average (centre weight 1/2)
                }
            }
        }

        // Rebuild h.x and derived surface arrays consistently (raise OR lower).
        for(int j = 0; j < jm; ++j){
            for(int k = 0; k < km; ++k){
                if(!massif[j][k]) continue;
                int ns = surf[j][k];
                if(ns < 1) ns = 1;                 // never clear an interior massif column to ocean
                const int old_s = i_topography[j][k];
                if(ns == old_s) continue;
                if(ns > old_s)
                    for(int i = old_s; i < ns; ++i) h.x[i][j][k] = 1.0;   // raise: add land
                else
                    for(int i = ns; i < old_s; ++i) h.x[i][j][k] = 0.0;   // lower: remove land
                i_topography[j][k] = ns;
                i_landscape[j][k]  = get_layer_height(ns);
                Landscape.y[j][k]  = get_layer_height(ns);
                Topography.y[j][k] = get_layer_height(ns);
            }
        }
    }



// rewriting bathymetrical data from -180° _ 0° _ +180° coordinate system to 0°- 360°
    #pragma omp parallel for schedule(dynamic)
    for(int j = 0; j < jm; j++){
        move_data(Topography.y[j], km);
        move_data(i_topography[j], km);
        move_data(i_landscape[j], km);
        move_data(Landscape.y[j], km);

        for(int i = 0; i < im; i++){
            move_data(h.x[i][j], km);
        }
    }


    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){

            int i_mount = i_topography[j][k];

            for(int i = 0; i < im; i++){

                if(i == i_mount)
                    Landscape.y[j][k] = get_layer_height(i);
                break;
            }
        }
    }

    cout << "      AGCM: init_topography ended" << endl;
}

// ==================== FULL-3D-STATE CHECKPOINT / RESTART ====================
// Binary dump of the prognostic fields so a debug run can resume at a chosen iter
// instead of re-spinning the dry circulation from scratch. Only the genuinely
// prognostic arrays are stored; densities, forces and all moist-convection
// diagnostics (S_*, M_*, MC_*, q_*_u, …) are recomputed each iteration from these.
std::vector<Array*> cAtmosphereModel::restart_arrays(){
    return { &t, &u, &v, &w, &c, &cloud, &ice, &gr, &co2,
             &tn, &un, &vn, &wn, &cn, &cloudn, &icen, &grn, &co2n,
             &p_dyn, &p_stat,
             &tke, &dis, &tken, &disn, &nue,
             &P_rain, &P_snow, &P_rainn, &P_snown };
}

void cAtmosphereModel::save_state(int iter, int Ma){
    const string fn = output_path + "/atm_restart_" + std::to_string(Ma) + "Ma_"
                      + std::to_string(iter) + ".bin";
    std::ofstream f(fn, std::ios::binary);
    if(!f){
        cout << "      AGCM: save_state FAILED to open " << fn << endl;
        return;
    }
    // Header: magic, dimensions, and the iter to resume at.
    const int32_t hdr[5] = { 0x41544D31 /*"ATM1"*/, im, jm, km, total_iter_count };
    f.write(reinterpret_cast<const char*>(hdr), sizeof(hdr));

    std::vector<Array*> arrs = restart_arrays();
    for(Array* a : arrs)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++)
                f.write(reinterpret_cast<const char*>(a->x[i][j]), km * sizeof(double));

    cout << "      AGCM: save_state wrote " << arrs.size() << " arrays to "
         << fn << " (total_iter_count " << total_iter_count << ")" << endl;
}

bool cAtmosphereModel::load_state(int iter, int Ma){
    const string fn = output_path + "/atm_restart_" + std::to_string(Ma) + "Ma_"
                      + std::to_string(iter) + ".bin";
    std::ifstream f(fn, std::ios::binary);
    if(!f){
        cout << "      AGCM: load_state: no file " << fn
             << " — running from scratch" << endl;
        return false;
    }
    int32_t hdr[5];
    f.read(reinterpret_cast<char*>(hdr), sizeof(hdr));
    if(!f || hdr[0] != 0x41544D31 || hdr[1] != im || hdr[2] != jm || hdr[3] != km){
        cout << "      AGCM: load_state: bad header / grid mismatch in " << fn
             << " — running from scratch" << endl;
        return false;
    }

    std::vector<Array*> arrs = restart_arrays();
    for(Array* a : arrs)
        for(int i = 0; i < im; i++)
            for(int j = 0; j < jm; j++){
                f.read(reinterpret_cast<char*>(a->x[i][j]), km * sizeof(double));
                if(!f){
                    cout << "      AGCM: load_state: truncated file " << fn
                         << " — running from scratch" << endl;
                    return false;
                }
            }

    total_iter_count = hdr[4];
    cout << "      AGCM: load_state restored " << arrs.size() << " arrays from "
         << fn << " (resuming at total_iter_count " << total_iter_count << ")" << endl;
    return true;
}
