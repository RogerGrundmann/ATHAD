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


// The Hadean surface is featureless — there is no topography to load. The variable
// keeps its inherited name because it is only an output-file stem downstream
// (UtilsAtm::writeFile, the ParaView writers, AtmospherePlotData) and renaming it
// would break the symbol-for-symbol correspondence with ATOM_Precipitation.
    bathymetry_name = "Hadean";

    init_topography();




// ATHAD: there is no surface data to read.
//
//  The modern Earth run anchored its initial state on observation and reconstruction:
//  the Scotese et al. (2021) global/equatorial/polar paleo-temperature curves, the NASA
//  surface temperature and precipitation fields, and — for Ma > 0 — an EarthByte/pygplates
//  reconstruction invoked as an external python script. None of that reaches the Hadean:
//  the Scotese curves stop at ~540 Ma and there is no observational field for 4.4 Ga.
//
//  ATHAD therefore PRESCRIBES its surface state instead of reading it. The surface
//  temperature comes from the t_surf_equator / t_surf_pole parameters, applied as a
//  pole-to-equator parabola in initTemperatureData(). See CLAUDE.md — this is the same
//  reason radiation cannot run in the inherited mode 5 (relaxation toward a Scotese target).

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
// ATHAD: read_Hydrosphere_SST() lived here — the atmosphere half of the atm<->ocean
// Picard loop, blending a previous hydrosphere run's sea-surface temperature into t.x[0].
// There is no hydrosphere in ATHAD and no ocean to couple to: the Hadean surface is a
// magma ocean whose temperature is prescribed, not solved by a companion model. Removed
// with its two parameters (sst_coupling_alpha, hyd_sst_iter).
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
void cAtmosphereModel::init_topography(){

    cout << endl << endl << endl << "      AGCM: init_topography" << endl;

//  ATHAD: the Hadean surface is featureless.
//
//  The topography of the Earth at ~4.4 Ga is unknown, so ATHAD prescribes a flat
//  global surface rather than guessing one. This is a physical choice, not a
//  missing input file: there is no bathymetry to read, no land/sea mask, no
//  coastline and no orography.
//
//  Everything the inherited Earth code derives from topography is pinned to its
//  ocean/sea-level default here:
//      h            = 0  -> AtomUtils::is_land() is false at every point, so all
//                           land branches in the RHS, BCs, turbulence, moist
//                           convection and ice schemes become dead by construction.
//      i_topography = 0  -> the surface is grid level 0 in every column.
//      Topography, Landscape, i_landscape = height of level 0 (i.e. 0 m).
//
//  Consequently the whole apparatus this function carried for the modern Earth —
//  the x-y-z file read, the peak/needle smoothing, the orographic slope cap and the
//  targeted massif smoothing, all of which existed to tame the Himalaya/Tibet
//  rampart that seeded a pressure-divergence runaway — has no subject and is gone.
//  Do not reintroduce it. See CLAUDE.md, invariant 1.

    h.initArray(im, jm, km, 0.0);

    const double surface_height = get_layer_height(0);                  // 0 m by construction of the stretched grid

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            i_topography[j][k] = 0;
            i_landscape[j][k]  = surface_height;
            Topography.y[j][k] = surface_height;
            Landscape.y[j][k]  = surface_height;
        }
    }

    cout << "      AGCM: flat Hadean surface — no topography, "
         << jm * km << " surface points, all water" << endl;

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
