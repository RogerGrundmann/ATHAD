// ============================================================================
// rad_selftest — standalone OFFLINE driver for MultiLayerRadiation.
//
// Exercises the REAL cAtmosphereModel::MultiLayerRadiation friend class on a
// single representative (US-standard) atmospheric column, WITHOUT wiring any
// call site into RunTimeSlice / run_3D_loop. Purpose: confirm the scheme runs
// on the model's actual 41-layer stretched grid with the p_stat units fix,
// produces no NaN, and yields physically sane epsilon / radiation / T + OLR.
//
// Build (from repo root, links the existing libatom.a — the friend line in
// cAtmosphereModel.h does not change the class layout so no rebuild is needed):
//   g++ -std=c++17 -march=native -fopenmp -Ilib -Iatmosphere -Ihydrosphere \
//       -Itinyxml2 test/rad_selftest.cpp -L. -latom -o rad_selftest
//   ./rad_selftest
// (compiled WITHOUT -ffast-math so the isnan() output check is reliable.)
// ============================================================================
#include "cAtmosphereModel.h"
#include "MultiLayerRadiation.h"
#include <cmath>
#include <cstdio>
#include <vector>

class RadiationSelfTest {
public:
    void run() {
        cAtmosphereModel m;                         // ctor: SetDefaultConfig + rad/the/phi alloc
        const int im = m.im, jm = m.jm, km = m.km;

        // radial coordinate, then layer heights. Was a third hardcoded copy of r0 = 1.0 and
        // dr = 0.025; it now takes whatever the model's single definition says, so the test keeps
        // exercising the real grid if that definition ever changes.
        m.initGridCoordinates();
        m.init_layer_heights();                     // private -> reachable via friend

        // allocate the fields the scheme reads / writes (ctor leaves these unallocated)
        m.t.initArray(im, jm, km, 1.0);
        m.c.initArray(im, jm, km, 0.0);
        m.co2.initArray(im, jm, km, 380.0);         // ppm (co2 is ppm-valued; 380 = 1x, 760 = 2x)
        m.p_stat.initArray(im, jm, km, 0.0);
        m.h.initArray(im, jm, km, 0.0);             // 0 everywhere -> ocean surface at i=0
        m.radiation.initArray(im, jm, km, 0.0);
        m.epsilon.initArray(im, jm, km, 0.0);
        m.cloud.initArray(im, jm, km, 0.0);         // cloud liquid mixing ratio [kg/kg]
        m.ice.initArray(im, jm, km, 0.0);           // cloud ice   mixing ratio [kg/kg]
        m.albedo.initArray_2D(jm, km, 0.0);
        m.epsilon_2D.initArray_2D(jm, km, 0.0);

        const int j0 = jm / 2, k0 = km / 2;           // equatorial-ish probe column

        // (Re)build a US-standard column (T, p, RH=0.5) over all (j,k) for a given CO2 factor,
        // run the REAL scheme, print the probe column + OLR, and return the surface T [K].
        // Called for 1x and 2x CO2 to check the surface energy balance now WARMS the surface.
        auto run_case = [&](double co2f, bool verbose, double qc_low, double qi_high) -> double {
            for (int i = 0; i < im; i++) {
                double zkm = m.get_layer_height(i) / 1000.0;
                double T   = (zkm < 11.0) ? 288.15 - 6.5 * zkm : 216.65;                 // K
                double p   = (zkm < 11.0) ? 1013.25 * pow(1.0 - 0.0065*zkm*1000.0/288.15, 5.256)
                                          : 226.32 * exp(-(zkm - 11.0)*1000.0/6341.6);    // hPa
                double esat = 6.1078 * exp(17.2694*(T-273.15)/(T-35.86));                 // hPa
                double cc   = m.ep * (0.5*esat) / p;                                      // kg/kg (RH=0.5)
                // optional test condensate: a low liquid cloud (1-2.5 km) and high cirrus ice (8-11 km)
                double qc = (qc_low  > 0.0 && zkm >= 1.0 && zkm <= 2.5) ? qc_low  : 0.0;  // kg/kg
                double qi = (qi_high > 0.0 && zkm >= 8.0 && zkm <= 11.0) ? qi_high : 0.0; // kg/kg
                for (int j = 0; j < jm; j++)
                    for (int k = 0; k < km; k++) {
                        m.t.x[i][j][k]      = T / m.t_0;   // nondim temperature (reset each case)
                        m.p_stat.x[i][j][k] = p;
                        m.c.x[i][j][k]      = cc;
                        m.co2.x[i][j][k]    = co2f;        // ppm: 380 -> 1x, 760 -> 2x
                        m.cloud.x[i][j][k]  = qc;
                        m.ice.x[i][j][k]    = qi;
                    }
            }
            std::vector<double> Tin(im);
            for (int i = 0; i < im; i++) Tin[i] = m.t.x[i][j0][k0] * m.t_0;

            MultiLayerRadiation(m).run();                  // the real scheme

            bool nan = false;
            if (verbose) printf("\n  i    z[m]   Tin[K]  Tout[K]   eps      rad[W/m2]\n");
            for (int i = 0; i < im; i++) {
                double Tout = m.t.x[i][j0][k0] * m.t_0, eps = m.epsilon.x[i][j0][k0], rad = m.radiation.x[i][j0][k0];
                if (std::isnan(Tout)||std::isnan(eps)||std::isnan(rad)||std::isinf(Tout)||std::isinf(eps)||std::isinf(rad)) nan = true;
                if (verbose && (i % 5 == 0 || i == im-1))
                    printf("%3d %7.0f  %6.1f  %7.2f  %.4f  %9.2f\n", i, m.get_layer_height(i), Tin[i], Tout, eps, rad);
            }
            double F = m.sigma * pow(m.t.x[0][j0][k0]*m.t_0, 4);
            for (int i = 1; i < im; i++) { double eps = m.epsilon.x[i][j0][k0]; F = F*(1.0-eps) + eps*m.sigma*pow(m.t.x[i][j0][k0]*m.t_0, 4); }
            double Tsurf = m.t.x[0][j0][k0]*m.t_0;
            printf("  CO2=%.0f ppm:  surface T = %.2f K (%.2f C)   OLR = %.1f W/m2   albedo = %.3f%s\n",
                   co2f, Tsurf, Tsurf-273.15, F, m.albedo.y[j0][k0], nan ? "   *** NaN/Inf ***" : "");
            return Tsurf;
        };

        double Ts1 = run_case(380.0, true, 0.0, 0.0);
        double Ts2 = run_case(760.0, false, 0.0, 0.0);
        printf("\n  ==> CO2 doubling surface warming: dT_s = %+.3f K   (should be POSITIVE)\n", Ts2 - Ts1);

        // Cloud/ice LW greenhouse: adding suspended condensate should LOWER OLR and WARM the
        // surface (clear -> low liquid cloud -> + high cirrus ice), at fixed CO2 = 380 ppm.
        printf("\n  --- cloud/ice LW optical depth (CO2 = 380 ppm) ---\n");
        printf("  clear sky:\n");   double Tsc  = run_case(380.0, false, 0.0,   0.0);
        printf("  low liquid cloud (qc=1e-4 kg/kg, 1-2.5 km):\n"); double TsL = run_case(380.0, false, 1.0e-4, 0.0);
        printf("  + high cirrus ice (qi=2e-5 kg/kg, 8-11 km):\n");  double TsI = run_case(380.0, false, 1.0e-4, 2.0e-5);
        printf("  ==> low-cloud surface LW warming: dT_s = %+.3f K ; +cirrus adds %+.3f K (both should be POSITIVE)\n",
               TsL - Tsc, TsI - TsL);

        // Cloud-CO2 overlap: does a thick cloud MASK the CO2 forcing? Compare the CO2 doubling
        // surface warming clear-sky vs with a thick low cloud present (grey-opaque clouds can
        // saturate the LW and shrink the CO2 marginal effect -> reduced sensitivity under cloud).
        printf("\n  --- cloud-CO2 overlap (CO2 doubling with vs without cloud) ---\n");
        double d_clear = run_case(760.0, false, 0.0, 0.0)    - run_case(380.0, false, 0.0, 0.0);
        double d_cloud = run_case(760.0, false, 6.9e-4, 0.0) - run_case(380.0, false, 6.9e-4, 0.0); // ~IC tropical LWC
        printf("  ==> CO2 doubling dT_s: clear = %+.3f K ; under thick low cloud = %+.3f K  (masking = %.0f%%)\n",
               d_clear, d_cloud, 100.0*(1.0 - d_cloud/d_clear));
    }
};

int main() { RadiationSelfTest().run(); return 0; }
