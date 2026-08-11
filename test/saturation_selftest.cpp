// Standalone self-test for SaturationH2O.h — checks the IAPWS correlations against
// published reference points before the model is allowed to depend on them.
//
// Build and run:
//   g++ -std=c++17 -I../atmosphere -I../lib -o saturation_selftest saturation_selftest.cpp && ./saturation_selftest
//
// or from the project root:  make test-saturation

#include "SaturationH2O.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

static int failures = 0;

static void check(const char* what, double got, double want, double rel_tol)
{
    const double rel = std::fabs(got - want) / std::fabs(want);
    const bool ok = (rel <= rel_tol);
    if (!ok) failures++;
    std::printf("  %-52s got %14.6g   want %14.6g   rel %8.2e   %s\n",
                what, got, want, rel, ok ? "ok" : "FAIL");
}

static void check_true(const char* what, bool cond)
{
    if (!cond) failures++;
    std::printf("  %-52s %s\n", what, cond ? "ok" : "FAIL");
}

int main()
{
    using namespace SaturationH2O;

    std::printf("\nSaturationH2O self-test\n\n");

    // ---- Saturation pressure over liquid water (IAPWS-95 reference points) ----
    std::printf(" saturation pressure over water [hPa]\n");
    // Triple point: 611.657 Pa exactly, by definition of the correlation's anchor.
    check("p_sat(273.16 K) = triple point",      saturationPressure(273.16),   6.11657,   2e-3);
    // Normal boiling point: 101.325 kPa at 373.124 K (the 1 atm boiling temperature).
    check("p_sat(373.124 K) = 1 atm",            saturationPressure(373.124),  1013.25,   2e-3);
    // Mid-range reference from the IAPWS-95 saturation table.
    check("p_sat(400 K)",                        saturationPressure(400.0),    2457.4,    5e-3);
    check("p_sat(500 K)",                        saturationPressure(500.0),   26392.0,    5e-3);
    check("p_sat(600 K)",                        saturationPressure(600.0),  123450.0,    5e-3);
    // Critical point.
    check("p_sat(647.096 K) = critical pressure",saturationPressure(647.095), 220640.0,   1e-3);

    // ---- Supercritical: the quantity must not exist ----
    std::printf("\n supercritical behaviour\n");
    check_true("p_sat(647.096 K) returns NO_SATURATION",
               saturationPressure(647.096) == NO_SATURATION);
    check_true("p_sat(1500 K) returns NO_SATURATION",
               saturationPressure(1500.0) == NO_SATURATION);
    check_true("isSupercritical(1500 K)", isSupercritical(1500.0));
    check_true("!isSupercritical(600 K)", !isSupercritical(600.0));

    // ---- Sublimation over ice (IAPWS 2011) ----
    std::printf("\n sublimation pressure over ice [hPa]\n");
    check("p_subl(273.16 K) = triple point",     sublimationPressure(273.16),  6.11657,   1e-4);
    check("p_subl(250 K)",                       sublimationPressure(250.0),   0.76053,   5e-3);
    check("p_subl(200 K)",                       sublimationPressure(200.0),   1.6362e-3, 1e-2);

    // ---- Monotonicity across the whole ATHAD range ----
    std::printf("\n monotonicity 150 K -> 647 K\n");
    bool mono = true;
    double prev = -1.0;
    for (double T = 150.0; T < AtmMixture::T_CRIT_H2O; T += 0.5) {
        const double p = saturationPressureAuto(T);
        if (!(p > prev) || !std::isfinite(p)) { mono = false; break; }
        prev = p;
    }
    check_true("saturationPressureAuto is finite and increasing", mono);

    // ---- Latent heat ----
    std::printf("\n latent heat of vaporisation [J/kg]\n");
    check("L(273.15 K)",                         latentHeat(273.15),          2.501e6,   1e-6);
    check("L(373.15 K) approx 2.257e6",          latentHeat(373.15),          2.257e6,   3e-2);
    check_true("L -> 0 at the critical point",   latentHeat(AtmMixture::T_CRIT_H2O) == 0.0);
    // Watson's 0.38 exponent takes L to zero slowly, so an absolute threshold near T_crit
    // is arbitrary. The meaningful statements are that L becomes a small FRACTION of its
    // cold value and that it is monotone all the way to zero.
    check_true("L < 2% of its 273 K value at T_crit - 0.006 K",
               latentHeat(647.09) < 0.02 * latentHeat(273.15));
    check_true("L decreases with T",             latentHeat(300.0) > latentHeat(500.0));
    check_true("L monotone decreasing to zero", [&]{
        double prev = latentHeat(200.0);
        for (double T = 200.5; T <= AtmMixture::T_CRIT_H2O; T += 0.5) {
            const double L = latentHeat(T);
            if (!(L <= prev) || L < 0.0 || !std::isfinite(L)) return false;
            prev = L;
        }
        return latentHeat(AtmMixture::T_CRIT_H2O) == 0.0;
    }());

    // ---- Exact vs dilute mass fraction ----
    // M_other: the ATHAD non-H2O mixture. CO2 and background at their configured ratio
    // give about 34 g/mol; the exact value does not matter for the structural checks.
    const double M_other = 0.034;
    std::printf("\n saturation mass fraction (exact conversion)\n");

    // Dilute limit: at Earth-like conditions the exact form must agree with the classic
    // dilute expression to well under a percent.
    {
        const double T = 288.0, p = 1013.25;
        const double E = saturationPressure(T);
        const double q_exact  = saturationMassFraction(E, p, 0.028964);
        const double ep       = 0.62198;                       // R_dry/R_vap for Earth air
        const double q_dilute = ep * E / (p - (1.0 - ep) * E);
        check("dilute limit: exact vs Magnus-style form", q_exact, q_dilute, 5e-3);
    }

    // Bulk limit: when the saturation pressure approaches the total pressure the mass
    // fraction must approach 1 smoothly, never exceed it, and never go negative.
    check("q_sat -> 1 as E -> p",   saturationMassFraction(999.0, 1000.0, M_other), 0.998114, 1e-4);
    check_true("q_sat = 1 when E > p",  saturationMassFraction(2000.0, 1000.0, M_other) == 1.0);
    check_true("q_sat = 1 when supercritical",
               saturationMassFraction(NO_SATURATION, 1000.0, M_other) == 1.0);
    check_true("q_sat in [0,1] over a wide sweep", [&]{
        for (double T = 150.0; T < 650.0; T += 1.0)
            for (double p = 1.0; p < 3.0e5; p *= 2.0) {
                const double q = saturationMassFractionAt(T, p, M_other);
                if (!(q >= 0.0 && q <= 1.0) || !std::isfinite(q)) return false;
            }
        return true;
    }());

    // ---- The condition that matters for ATHAD ----
    std::printf("\n ATHAD condensation level\n");
    // At the surface: 1500 K, 250 bar -> supercritical, nothing condenses.
    check_true("surface (1500 K, 250 bar): no condensation limit",
               saturationMassFractionAt(1500.0, 250000.0, M_other) == 1.0);
    // A point where saturation genuinely bites. Note 640 K / 20 bar does NOT: the
    // saturation pressure there is 186 bar, so water at 20 bar is superheated and q_sat
    // is correctly 1 (no constraint). Condensation needs p_H2O to exceed p_sat(T), which
    // takes both a lower temperature and a high total pressure.
    {
        const double q = saturationMassFractionAt(500.0, 100000.0, M_other);
        std::printf("  q_sat(500 K, 100 bar) = %.4f  (Hadean q_H2O = 0.6724)\n", q);
        check_true("q_sat < Hadean q_H2O at 500 K / 100 bar", q < 0.6724);
        check_true("q_sat > 0",                              q > 0.0);
    }
    {
        const double q = saturationMassFractionAt(640.0, 20000.0, M_other);
        std::printf("  q_sat(640 K, 20 bar)  = %.4f  (superheated: p_sat = 186 bar)\n", q);
        check_true("superheated vapour carries no condensation limit", q == 1.0);
    }

    std::printf("\n%s  (%d failure%s)\n\n",
                failures ? "SELF-TEST FAILED" : "self-test passed",
                failures, failures == 1 ? "" : "s");
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
