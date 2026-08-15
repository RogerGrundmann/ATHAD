#pragma once

#include "cAtmosphereModel.h"
#include "Utils.h"

#include <iostream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace AtomUtils;

class VelocityInitializer {
public:
    explicit VelocityInitializer(cAtmosphereModel& model)
        : m(model)
    {}

    void compute()
    {
        using namespace std;
        cout << endl << "      AGCM: init_velocities" << endl;
        if(latScale() != 1.0)
            cout << "      cell latitudes scaled by " << latScale()
                 << ": Hadley anchor at " << (90 - js(75)) << " deg, Ferrel at "
                 << (90 - js(45)) << " deg, polar at " << (90 - js(15)) << " deg" << endl;

        // u-component up to tropopause and back on half distance
        init_u(m.u, js(0));
        init_u(m.u, js(30));
        init_u(m.u, js(60));
        init_u(m.u, js(90));
        init_u(m.u, js(120));
        init_u(m.u, js(150));
        init_u(m.u, js(180));

        // initialise v: tropopause and surface values per latitude
        // equator
        init_v_or_w(m.v, js(90),  0.0,  0.0);                              // lat:   0   j=90
        // northern polar cell
        init_v_or_w(m.v, js(0),  0.5,  0.0);                              // lat:  90   j=0
        init_v_or_w(m.v, js(15),  0.5,  0.6);                              // lat:  75   j=15
        // southern polar cell
        init_v_or_w(m.v, js(180),  0.5,  0.0);                              // lat: -90   j=180
        init_v_or_w(m.v, js(165),  0.5,  0.6);                              // lat: -75   j=165
        // northern Ferrel cell
        init_v_or_w(m.v, js(30), -0.2,  0.0);                              // lat:  60   j=30
        init_v_or_w(m.v, js(45),  4.0, -1.5);                              // lat:  45   j=45
        // southern Ferrel cell
        init_v_or_w(m.v, js(150), -0.2,  0.0);                              // lat: -60   j=150
        init_v_or_w(m.v, js(135),  4.0, -1.5);                              // lat: -45   j=135
        // Hadley cells — SYMMETRIC about the equator.
        //
        // These two carried 4.0 at 15N against 3.0 at 15S (and, before someone swapped
        // them, 3.0 against 4.0 — the commented-out pair that used to sit here). Every
        // other mirror pair in this routine is identical, so this was the only asymmetry
        // in the whole velocity initialisation, and it showed: the meridional
        // streamfunction came out with the southern cell stronger than the northern by
        // 366/278 = 1.32 at 15 degrees, against 4.0/3.0 = 1.33, and 11 % antisymmetry
        // error globally that had not washed out by iteration 100.
        //
        // On Earth a north-south Hadley asymmetry is physical: the ITCZ sits north of the
        // equator because of the land-sea distribution. ATHAD has no land, no topography,
        // a symmetric prescribed surface temperature, an insolation profile that is
        // EXPLICITLY mirrored (short_wave_radiation[j] = short_wave_radiation[j_max-j]),
        // no obliquity and no seasons. Nothing here can sustain a hemispheric asymmetry,
        // so all of it was inherited from these two numbers.
        //
        // 3.5 is their mean, which removes the asymmetry and leaves the total initial
        // Hadley mass flux unchanged.
        init_v_or_w(m.v, js(60),  0.0,  0.5);                              // lat:  30   j=60
        init_v_or_w(m.v, js(75), -3.0,  3.5);                              // lat:  15   j=75
        init_v_or_w(m.v, js(120),  0.0,  0.5);                              // lat: -30   j=120
        init_v_or_w(m.v, js(105), -3.0,  3.5);                              // lat: -15   j=105

        // initialise w: tropopause and surface values per latitude.
        // w is the ZONAL jet (East+). The SURFACE value (2nd coeff) is what the
        // atm->ocean transfer hands to the ocean, so it must reproduce the observed
        // surface wind BANDS: trade EASTERLIES (w<0) through the tropics/subtropics
        // (0-30deg, peak ~15deg), mid-latitude WESTERLIES (w>0) peaking ~45deg, then
        // weakening poleward. The wind-stress CURL between the trade easterlies and
        // the mid-latitude westerlies is what drives the subtropical (anticyclonic)
        // gyre; the curl between the westerly max and the pole drives the subpolar
        // (cyclonic) gyre. The tropopause value (1st coeff) keeps the upper-level
        // westerly jets (subtropical jet strongest at 30deg). See wind-IC diagnosis
        // in project_hydro_ekman_sh_gyre.
        // equator
        init_v_or_w(m.w, js(90), -3.0, -5.0);                             // lat:   0   j=90   easterly (equatorial)
        // northern polar cell
        init_v_or_w(m.w, js(0),  0.0,  0.0);                              // lat:  90   j=0
        // southern polar cell
        init_v_or_w(m.w, js(180),  0.0,  0.0);                              // lat: -90   j=180
        // northern Ferrel cell (mid-latitude westerlies, weakening to the pole)
        init_v_or_w(m.w, js(30), 10.0,  6.0);                             // lat:  60   j=30   westerly
        // southern Ferrel cell
        init_v_or_w(m.w, js(150), 10.0,  6.0);                             // lat: -60   j=150  westerly
        // northern subtropics — horse latitudes (trade/westerly transition, calm)
        init_v_or_w(m.w, js(60), 30.0, -1.0);                             // lat:  30   j=60   weak easterly
        // southern subtropics
        init_v_or_w(m.w, js(120), 30.0, -1.0);                             // lat: -30   j=120  weak easterly
        // northern westerly max at j=45
        init_v_or_w(m.w, js(45), 15.0, 10.0);                             // lat:  45   j=45   westerly max
        // southern westerly max at j=135
        init_v_or_w(m.w, js(135), 15.0, 10.0);                             // lat: -45   j=135  westerly max
        // northern trade-easterly max at j=75 (15N)
        init_v_or_w(m.w, js(75),  5.0, -7.0);                             // lat:  15   j=75   easterly (trade max)
        // southern trade-easterly max at j=105 (15S)
        init_v_or_w(m.w, js(105),  5.0, -7.0);                             // lat: -15   j=105  easterly (trade max)

        // forming diagonals — northern hemisphere
        form_diagonals(m.u, js(0), js(30));
        form_diagonals(m.w, js(0), js(30));
        form_diagonals(m.w, js(30), js(45));
        form_diagonals(m.v, js(0), js(15));
        form_diagonals(m.v, js(15), js(30));

        form_diagonals(m.u, js(30), js(60));
        form_diagonals(m.w, js(45), js(60));
        form_diagonals(m.v, js(30), js(45));
        form_diagonals(m.v, js(45), js(60));

        form_diagonals(m.u, js(60), js(90));
        form_diagonals(m.w, js(60), js(75));                                   // 30N->15N (trade node at j=75)
        form_diagonals(m.w, js(75), js(90));                                   // 15N->0
        form_diagonals(m.v, js(60), js(75));
        form_diagonals(m.v, js(75), js(90));

        // forming diagonals — southern hemisphere
        form_diagonals(m.u, js(90), js(120));
        form_diagonals(m.w, js(90), js(105));                                  // 0->15S (trade node at j=105)
        form_diagonals(m.w, js(105), js(120));                                 // 15S->30S
        form_diagonals(m.w, js(120), js(135));
        form_diagonals(m.v, js(90), js(105));
        form_diagonals(m.v, js(105), js(120));

        form_diagonals(m.u, js(120), js(150));
        form_diagonals(m.w, js(135), js(150));
        form_diagonals(m.v, js(120), js(135));
        form_diagonals(m.v, js(135), js(150));

        form_diagonals(m.u, js(150), js(180));
        form_diagonals(m.w, js(150), js(180));
        form_diagonals(m.v, js(150), js(165));
        form_diagonals(m.v, js(165), js(180));

        // Zero land cells; non-dimensionalise air cells — single fused pass
        const double inv_u_0 = 1.0 / m.u_0;

        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < m.im; i++) {
            for (int k = 0; k < m.km; k++) {
                for (int j = 0; j < m.jm; j++) {
                    if (is_land(m.h, i, j, k)) {
                        m.u.x[i][j][k] = 0.0;
                        m.v.x[i][j][k] = 0.0;
                        m.w.x[i][j][k] = 0.0;
                    } else {
                        m.u.x[i][j][k] *= inv_u_0;
                        m.w.x[i][j][k] *= inv_u_0;
                        if (!m.use_NASA_velocity && j > 90) {
                            m.v.x[i][j][k] = -m.v.x[i][j][k] * inv_u_0;
                        } else {
                            m.v.x[i][j][k] *= inv_u_0;
                        }
                    }
                }
            }
        }
/*
        // Surface taper on v and w (the two HORIZONTAL components in this model's
        // (r,θ,φ) convention: v = meridional, w = zonal): linearly damp from the
        // local value at i=5 down to zero at i=0, so the lowest five layers carry no
        // horizontal wind at the ground reference and grow smoothly into the
        // prescribed profile above. u is left alone because u is the RADIAL/VERTICAL
        // velocity here (NOT the zonal jet — that is w); init_u already gives it a
        // small profile that ramps to zero at the surface.
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < m.jm; j++) {
            for (int k = 0; k < m.km; k++) {
                for (int i = 0; i <= 5; i++) {
                    const double factor = static_cast<double>(i) / 5.0;
                    m.v.x[i][j][k] *= factor;
                    m.w.x[i][j][k] *= factor;
                }
            }
        }
*/
    cout << "      AGCM: init_velocities ended" << endl;
    }

    // Linear blend of u/v/w across j in [lat-3, lat+3].
    // Currently not called by compute() (dead code in the original),
    // but kept here as it logically belongs with velocity initialisation.
    void smooth_transition(int lat)
    {
        const int    start     = lat - 3;
        const int    end       = lat + 3;
        const double inv_range = 1.0 / (double)(end - start);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                const double u_start = m.u.x[i][start][k];
                const double v_start = m.v.x[i][start][k];
                const double w_start = m.w.x[i][start][k];
                const double u_slope = (m.u.x[i][end][k] - u_start) * inv_range;
                const double v_slope = (m.v.x[i][end][k] - v_start) * inv_range;
                const double w_slope = (m.w.x[i][end][k] - w_start) * inv_range;
                for (int j = start; j <= end; j++) {
                    const double t     = (double)(j - start);
                    m.u.x[i][j][k] = u_slope * t + u_start;
                    m.v.x[i][j][k] = v_slope * t + v_start;
                    m.w.x[i][j][k] = w_slope * t + w_start;
                }
            }
        }
    }

private:
    cAtmosphereModel& m;


    // ATM_CELL_LAT_SCALE — compress the prescribed cell latitudes toward the equator.
    //
    // Every anchor below sits at an EARTH latitude: Hadley 15 deg, Ferrel 45, polar 75,
    // with the trade/westerly nodes between them. Those latitudes are a consequence of
    // Earth's thermal Rossby number, and this atmosphere's is 12.5x smaller —
    // Ro_T = g*H*(dtheta/theta)/(omega^2 a^2) = 0.0048 against 0.0598, because rotation is
    // 4.35x faster (omega^2 18.9x) and the FRACTIONAL equator-pole contrast is 4.7x weaker
    // (50 K on 1500 K against 45 K on 288 K), only partly offset by a 7x deeper atmosphere.
    // Held-Hou then puts the direct cell's edge near 5 deg here against 18 deg for Earth,
    // i.e. cells roughly a third as wide. A hot surface is not a strongly DIFFERENTIALLY
    // heated one, which is why the higher energy content narrows the circulation instead of
    // widening it.
    //
    // The knob exists to test whether the cell decay measured in README item 28 is partly
    // the model rejecting an over-wide imposed cell, rather than only the residual
    // meridional force. 1.0 = Earth's latitudes (default, what every run so far used).
    // The poles are fixed points of the map, so the interpolation still spans [0, jm-1].
    static double latScale(){
        static const double s = [](){ const char* e = getenv("ATM_CELL_LAT_SCALE");
            return e ? atof(e) : 1.0; }();
        return s;
    }
    int js(int j) const {
        if(j <= 0) return 0;
        if(j >= m.jm - 1) return m.jm - 1;
        const int jeq = (m.jm - 1) / 2;
        return (int)std::lround(jeq + (double)(j - jeq) * latScale());
    }

    void init_u(Array& u, int j)
    {
        const double ua_00  = 0.02894;
        const double ua_30  = 0.02315;
        const double ua_60  = 0.01736;
        const double ua_90  = 0.011574;

        double coeff;
        switch (j) {
            case  90: coeff =  ua_00; break;
            case  60: coeff = -ua_30; break;
            case 120: coeff = -ua_30; break;
            case  30: coeff =  ua_60; break;
            case 150: coeff =  ua_60; break;
            case   0: coeff = -ua_90; break;
            case 180: coeff = -ua_90; break;
            default:  return;
        }

        const int    tl           = m.get_tropopause_layer(j);
        const double tropo_h      = m.get_layer_height(tl);
        const double half_tropo_h = tropo_h / 3.0;
        const double inv_ascent   = 3.0 / half_tropo_h;
        const double inv_descent  = 1.0 / half_tropo_h;

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < tl; i++) {
                const double h     = m.get_layer_height(i);
                const double ratio = (h < half_tropo_h)
                    ? h * inv_ascent
                    : (tropo_h - h) * inv_descent;
                u.x[i][j][k] = coeff * ratio;
            }
        }
    }

    void init_v_or_w(Array& v_or_w, int j, double coeff_trop, double coeff_sl)
    {
        const int    tl          = m.get_tropopause_layer(j);
        const double inv_tropo_h = 1.0 / m.get_layer_height(tl);

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            double sl = coeff_sl;
            if (m.use_NASA_velocity && is_ocean_surface(m.h, 0, j, k)) {
                sl = v_or_w.x[0][j][k];
            }
            const double slope = (coeff_trop - sl) * inv_tropo_h;

            for (int i = 0; i < tl; i++) {
                v_or_w.x[i][j][k] = slope * m.get_layer_height(i) + sl;
            }
        }

        init_v_or_w_above_tropopause(v_or_w, j, coeff_trop);
    }

    void init_v_or_w_above_tropopause(Array& v_or_w, int j, double coeff)
    {
        const int tl = m.get_tropopause_layer(j);
        if (tl >= m.im - 1) return;

        const double h_top     = m.get_layer_height(m.im - 1);
        const double inv_range = 1.0 / (h_top - m.get_layer_height(tl));

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = tl; i < m.im; i++) {
                v_or_w.x[i][j][k] = coeff * (h_top - m.get_layer_height(i)) * inv_range;
            }
        }
    }

    void form_diagonals(Array& a, int start, int end)
    {
        const double inv_range = 1.0 / (double)(end - start);

        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < m.km; k++) {
            for (int i = 0; i < m.im; i++) {
                const double a_start = a.x[i][start][k];
                const double slope   = (a.x[i][end][k] - a_start) * inv_range;
                for (int j = start; j < end; j++) {
                    a.x[i][j][k] = slope * (double)(j - start) + a_start;
                }
            }
        }
    }
};
