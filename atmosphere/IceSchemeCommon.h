#pragma once

#include "cAtmosphereModel.h"

// ============================================================================
// Shared machinery for the ice / precipitation microphysics schemes
// (Zero/One/Two/ThreeCatIceScheme). The four schemes share a large amount of
// physics-independent boilerplate — boundary extrapolation, topography fill,
// saturation, the finiteness flux cap — which was previously duplicated (and
// independently mis-fixed) in each file. It is factored out here so a fix lives
// in ONE place. The per-scheme MICROPHYSICS stays in each scheme's
// computeColumns(); only the common scaffolding is shared.
// ============================================================================
namespace IceSchemeCommon {

    // Hard cap on a precipitation flux [kg/(m2*s)] (~260 mm/d, well above any
    // physical precip). Backstop against the ∝Snow riming/deposition runaways.
    constexpr double P_max_flux = 3.0e-3;

    // c43/c13 quadratic edge extrapolation of ONE precipitation-flux field over
    // the top (radial), latitude (theta) and longitude (phi) boundaries. The phi
    // seam optionally gets periodicity averaging: rain uses it in every scheme;
    // snow/graupel use it in the >=2-category schemes but NOT in ZeroCat, so it
    // is a per-field flag to preserve each scheme's exact behaviour.
    inline void extrapolateBC(cAtmosphereModel& m, Array& P, bool periodic_avg) {
        // Top boundary
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++)
            for(int k = 0; k < m.km; k++)
                P.x[m.im-1][j][k] = m.c43 * P.x[m.im-2][j][k] - m.c13 * P.x[m.im-3][j][k];

        // Latitude (theta) boundaries
        #pragma omp parallel for collapse(2)
        for(int k = 0; k < m.km; k++)
            for(int i = 0; i < m.im; i++){
                P.x[i][0][k]      = m.c43 * P.x[i][1][k]      - m.c13 * P.x[i][2][k];
                P.x[i][m.jm-1][k] = m.c43 * P.x[i][m.jm-2][k] - m.c13 * P.x[i][m.jm-3][k];
            }

        // Longitude (phi) boundaries — von Neumann + optional periodicity average
        #pragma omp parallel for collapse(2)
        for(int i = 0; i < m.im; i++)
            for(int j = 0; j < m.jm; j++){
                P.x[i][j][0]      = m.c43 * P.x[i][j][1]      - m.c13 * P.x[i][j][2];
                P.x[i][j][m.km-1] = m.c43 * P.x[i][j][m.km-2] - m.c13 * P.x[i][j][m.km-3];
                if(periodic_avg)
                    P.x[i][j][0] = P.x[i][j][m.km-1] =
                        (P.x[i][j][0] + P.x[i][j][m.km-1]) * 0.5;
            }
    }

    // Fill the sub-surface land cells (i below i_topography) of ONE field with the
    // surface (i_mount) value. Each field independently copies its own ground value
    // downward, so calling this per field is equivalent to the old snapshot-all form.
    inline void fillTopography(cAtmosphereModel& m, Array& F) {
        #pragma omp parallel for collapse(2)
        for(int j = 0; j < m.jm; j++)
            for(int k = 0; k < m.km; k++){
                int i_mount = m.i_topography[j][k];
                double v = F.x[i_mount][j][k];
                for(int i = i_mount - 1; i >= 0; i--)
                    if(AtomUtils::is_land(m.h, i, j, k)) F.x[i][j][k] = v;
            }
    }

    // Saturation specific humidity over water / ice [kg/kg] (Magnus form used by
    // every scheme). Available for computeColumns to use; kept identical to the
    // inline expressions so a scheme can adopt it without changing results.
    inline double qSatWater(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        const double E = m.hp * AtomUtils::exp_func(t_u, 17.2694, 35.86);
        return m.ep * E / (m.p_stat.x[i][j][k] - E);
    }
    inline double qSatIce(cAtmosphereModel& m, double t_u, int i, int j, int k) {
        const double E = m.hp * AtomUtils::exp_func(t_u, 21.8746, 7.66);
        return m.ep * E / (m.p_stat.x[i][j][k] - E);
    }

    // ---- Vapour -> ice -> snow throttle (Seifert-Beheng / COSMO; TwoCat's) ----
    // The stable way to make snow: grow it through the CLOUD-ICE reservoir, not
    // directly from vapour. Deposition onto ice is SUPERSATURATION-LIMITED
    // (S_i_dep proportional to c - q_Ice, and to the ice particle population
    // N_i*m_i^(1/3)), and the ice->snow autoconversions (S_i_au aggregation,
    // S_d_au deposition) are bounded — so snow cannot run away. This is exactly
    // what OneCat lacked (it deposited vapour straight onto snow proportional to
    // Snow^0.58, an unbounded feedback). TwoCat carries the same physics inline;
    // sharing it here lets OneCat inherit the throttle.
    // Also returns the ice particle number density N_i and mean mass m_i, which the
    // callers need for their own nucleation / ice-rain-collection terms.
    struct IceSnowRates { double S_i_dep = 0.0, S_i_au = 0.0, S_d_au = 0.0,
                                 N_i = 0.0, m_i = 1.0e-9; };

    inline IceSnowRates depositionThrottle(cAtmosphereModel& m, int i, int j, int k,
                                           double t_u, double q_Ice, double dt_snow_dim) {
        constexpr double N_i_0   = 1.0e2;    // 1/m3
        constexpr double m_i_0   = 1.0e-12;  // kg
        constexpr double m_i_max = 1.0e-9;   // kg
        constexpr double m_s_0   = 3.0e-9;   // kg (snow-size threshold)
        constexpr double c_i_dep = 1.3e-5;   // m3/(s*kg^1/3)
        constexpr double c_i_au  = 1.0e-3;   // 1/s
        constexpr double t_hn    = 236.15;   // K (-37C) homogeneous freezing floor

        IceSnowRates r;
        const double c   = m.c.x[i][j][k];
        const double ice = m.ice.x[i][j][k];

        // ice number density + mean particle mass (only in the mixed-phase band)
        double N_i = 0.0, m_i = m_i_max;
        if(t_u <= m.t_0 && t_u > t_hn){
            N_i = N_i_0 * std::exp(0.2 * (m.t_0 - t_u));
            m_i = std::min(m.r_humid.x[i][j][k] * ice / N_i, m_i_max);
            m_i = std::max(m_i_0, std::min(m_i, m_i_max));
        }

        // deposition (supersaturation-limited) / sublimation (ice-mass-limited)
        if(c > q_Ice)      r.S_i_dep = c_i_dep * N_i * std::pow(m_i, 1.0/3.0) * (c - q_Ice);
        else if(c < q_Ice) r.S_i_dep = std::max(-ice / dt_snow_dim, (c - q_Ice) / dt_snow_dim);

        // ice -> snow: aggregation (bounded) + depositional autoconversion
        if(ice > 0.0)       r.S_i_au = std::max(c_i_au * ice, 0.0);
        if(r.S_i_dep > 0.0) r.S_d_au = r.S_i_dep / (1.5 * (std::pow(m_s_0 / m_i, 2.0/3.0) - 1.0));
        r.N_i = N_i;  r.m_i = m_i;
        return r;
    }
}
