# Given a parameter definition, generates necessary C++, Python and XML bindings
# coding=utf-8



def main():

    # read the input definition
    # name, description, datatype, default in the tuples


    PARAMS = {                                                          # dictionary{} (PARAMS) with keys ('common', etc.):[ and their tuples ('output_path', etc.)]
        'common': [
            # Relative to the directory the run is launched from. The ATOM line's convention
            # is to run from python/ and write to python/output_<name>/ — ATOM itself uses
            # output_ATOM/ — so a default of output_Hadean/ lands in python/output_Hadean.
            # (The giant-planet siblings hyphenate, output-Uranus/, but ATHAD is from the
            # ATOM line and follows it.)
            ('output_path', 'directory where model outputs should be placed(must end in /)', 'string', 'output_Hadean/'),



            ('config_xml_path', 'directory where the configuration files are located', 'string', '../python'),

            ('verbose', 'some description of the module', 'bool', False),

#            ('paraview_panorama_vts_flag','flag to control if create paraview panorama', 'bool', False),
            ('paraview_panorama_vts_flag','flag to control if create paraview panorama', 'bool', True),

            ('Coriolis', 'Coriolis force', 'double', 1),
#            ('Coriolis', 'Coriolis force', 'double', 0),

            ('centrifugal', 'centrifugal force', 'double', 1),
#            ('centrifugal', 'centrifugal force', 'double', 0),

            ('buoyancy', 'buoyancy force', 'double', 1),
#            ('buoyancy', 'buoyancy force', 'double', 0),

            ('r_Earth', 'radius of the Earth in km', 'double', 6370.001),
            # ATHAD: the Hadean day was far shorter — tidal braking has been slowing the
            # Earth ever since. A 5.5 h day at ~4.4 Ga gives omega = 2*pi/19800 s =
            # 3.17e-4 rad/s, 4.35x the modern 7.292e-5. ASSUMPTION: estimates for the
            # early rotation range over roughly 4-6 h, and the Rossby number scales with
            # it, so this directly sets how many circulation cells the model can support.
            ('omega', 'ATHAD: rotation rate in rad/s (5.5 h Hadean day)', 'double', 3.17e-4),

            ('g', 'gravitational acceleration of the earth in m/s²', 'double', 9.8066),








            # ATHAD is ONE EPOCH. The paleo-reconstruction inputs that lived here — the
            # topography grids, the NASA surface fields, the Scotese temperature curves, the
            # pygplates reconstruction script, the Ma switches and the hydrosphere SST
            # coupling — are all gone: nothing reconstructs 4.4 Ga, and there is no
            # hydrosphere to couple to. time_start/end/step remain only because the
            # time-slice loop is still structural; a single slice is all that runs.

            ('time_start', 'start time', 'int', 0),
#            ('time_end', 'end time', 'int', 10),
            ('time_end', 'end time', 'int', 0),
            ('time_step', 'step size between timeslices', 'int', 10),







            ('use_earthbyte_reconstruction', 'control whether use earthbyte method to recontruct grids', 'bool', False),
#            ('use_earthbyte_reconstruction', 'control whether use earthbyte method to recontruct grids', 'bool', True),

            ('use_NASA_velocity', 'if use NASA velocity to initialise velocity', 'bool', False),
#            ('use_NASA_velocity', 'if use NASA velocity to initialise velocity', 'bool', True),

            ('use_NASA_temperature', 'if use NASA temperature to initialise surface temperature', 'bool', True),



#            ('CategoryIceScheme', 'number chooses Three(3)-Category Ice Scheme with rain, snow and graupel', 'int', 3),
            ('CategoryIceScheme', 'number chooses Two(2)-Category Ice Scheme with rain, snow', 'int', 2),
#            ('CategoryIceScheme', 'number chooses One(1)-Category Ice Scheme with rain, snow', 'int', 1),
#            ('CategoryIceScheme', 'number chooses Zero(0)-Category Ice Scheme with rain (Warm Rain Scheme)', 'int', 0),
#            ('CategoryIceScheme', 'number chooses no scheme(-1) no precipitation', 'int', -1),

            # ATHAD: 250 bar = 250000 hPa. NOTE the surface pressure is not actually taken
            # from p_0 — InitValues/ThermoAtm build it as 1e-2*(r_air*R_Air*T), i.e. from the
            # reference density — so r_air below must be consistent with this value.
            ('p_0', 'pressure at sea level in hPa', 'double', 250000.0),
            ('t_0', 'temperature in K compare to 0°C', 'double', 273.15),

            # ATHAD: the Hadean surface temperature is PRESCRIBED, not read from a
            # paleo-temperature curve — no Scotese reconstruction reaches 4.4 Ga. A
            # runaway steam atmosphere at 250 bar is optically thick enough that the
            # equator-pole contrast is small; 50 K is an assumption, not a result.
            ('t_surf_equator', 'ATHAD: prescribed Hadean surface temperature at the equator in K', 'double', 1500.0),
            ('t_surf_pole', 'ATHAD: prescribed Hadean surface temperature at the poles in K', 'double', 1450.0),

            # ATHAD: the COSMO barometric profile T(h) = T0*sqrt(1 - 2*beta*g*h/(R*T0^2)) has a
            # near-surface lapse rate beta*g/(R*T0), so beta is derived rather than fixed:
            #
            #     beta = cosmo_lapse_fraction * R_mix * T_surf / cp
            #
            # cosmo_lapse_fraction is the lapse as a fraction of the DRY ADIABAT g/cp.
            #
            # Earth's value is 0.51: its observed ~5 K/km sits about half way to the 9.8 K/km
            # dry adiabat. That deficit is not a coincidence or a tuning constant — it is
            # LATENT HEAT. A rising saturated parcel condenses water and releases heat, which
            # partly offsets adiabatic cooling and flattens the lapse toward the moist adiabat.
            #
            # ATHAD's deep column condenses NOTHING: it is supercritical from the surface to
            # ~177 km, so no latent heat is released and there is nothing to flatten the lapse.
            # The physically consistent value is therefore 1.0, the dry adiabat of the mixture,
            # g/cp = 4.81 K/km. Carrying Earth's 0.51 across gave 2.46 K/km, a column that only
            # reached 1279 K at 300 km, never crossed the critical temperature, and so could
            # never condense — a self-fulfilling assumption, since the flattened lapse was
            # itself justified by condensation that the flattening then prevented.
            #
            # On the dry adiabat the column crosses 647 K at 177 km and reaches saturation
            # (p_H2O > p_sat) near 300 K at ~249 km and 0.04 bar, which is where a runaway
            # greenhouse should put its cloud deck.
            ('cosmo_lapse_fraction', 'ATHAD: near-surface lapse rate as a fraction of the dry adiabat; 1.0 = dry adiabat', 'double', 1.0),

            # ATHAD: temperature of the optically thin upper atmosphere, which is radiative
            # rather than convective. Set from the ENERGY BUDGET, not from a measured OLR:
            #
            #     sigma*T_skin^4 = (1 - albedo)*SW + geothermal
            #                    = 0.92*93.5 + 150 = 236.0 W/m2  ->  254.0 K
            #
            # An earlier value of 231.3 K came from T_skin = (OLR/2sigma)^(1/4) using a
            # previously MEASURED OLR, which is circular — the model then reproduced an OLR of
            # sigma*T_skin^4 and confirmed nothing but its own arithmetic.
            #
            # This is now the STARTING value of an iterated fixed point, not a fixed input.
            # With t_skin_relax > 0 the model re-derives it in the loop against its OWN mean
            # albedo — the clear-sky value here is only where the iteration begins.
            #
            # READ THIS BEFORE QUOTING THE OLR. The column's top is isothermal AT t_skin by
            # construction (ThermoAtm::densities re-imposes the profile every iteration) and
            # it is optically thick, so the emission level sits there and the reported OLR is
            # identically sigma*t_skin^4. Measured: t_skin = 254 K -> OLR 236.01 W/m2;
            # t_skin = 240 K -> OLR 188.13 W/m2; sigma*T^4 = 236.01 and 188.13. The OLR is an
            # INPUT wearing an output's clothes. Closing the fixed point below makes it equal
            # the absorbed flux — which is then true by construction, not by test. Making the
            # OLR a genuine prediction requires letting the top find its own temperature
            # radiatively instead of having it prescribed. See the README.
            #
            # It replaces the inherited t_00 = 236.15 K ("-37 C"), which was Earth's
            # tropopause temperature and landed near the right range here by coincidence.
            ('t_skin', 'ATHAD: starting radiative-equilibrium temperature of the optically thin top, in K', 'double', 254.0),

            # ATHAD: under-relaxation of the t_skin fixed-point iteration, per radiation call.
            #
            #     sigma*t_skin^4 = (1 - albedo_mean)*SW_mean + geothermal_flux
            #
            # with albedo_mean the model's OWN cos(latitude)-weighted albedo, built by
            # MultiLayerRadiation from the condensate the model actually made. Set to 0 to
            # hold t_skin at the configured value (the old behaviour, which left the budget
            # open by -32 W/m2). 0.25 converges in a handful of iterations and is well inside
            # the stability limit, the map being a gentle T^(1/4).
            ('t_skin_relax', 'ATHAD: relaxation of the t_skin fixed point per radiation call; 0 = hold t_skin fixed', 'double', 0.25),
            # ATHAD reference density of the MIXTURE at the surface, not of dry air:
            # rho = p/(R_mix*T) = 25e6 Pa / (387.9 * 1500 K) = 42.97 kg/m³ (Earth: 1.2041).
            # This is what sets the surface pressure, via p = 1e-2*(r_air*R_Air*T).
            ('r_air', 'ATHAD: reference density of the atmospheric mixture at the surface in kg/m³', 'double', 42.97),
            ('r_0_water', 'reference density of fresh water in kg/m3', 'double', 997.0),
            ('t_equat_modern', 'mean temperature of the modern earth in °C', 'double', 15.4),

 
            # ATHAD insolation: the faint young Sun, at the TOP OF THE ATMOSPHERE.
            #
            # At 4.4 Ga the solar constant was about 0.71 of today's 1361 W/m2, i.e.
            # S = 966.3 W/m2, so the global mean incident flux is S/4 = 241.6 W/m2.
            #
            # The previous values, 116 and 71, were ATOM_Precipitation's 163.3 and 100.0
            # scaled by 0.71 — and those are Earth's absorbed-at-the-SURFACE shortwave,
            # already reduced by Earth's albedo and its atmosphere's absorption.
            # MultiLayerRadiation uses these as the flux INCIDENT at the top and then
            # applies ATHAD's own albedo, so Earth's albedo was being counted twice: the
            # cos(latitude)-weighted mean came to 107.5 W/m2, lighting the planet at 45 %
            # of its own insolation.
            #
            # The values below are a fit of the model's parabola short_wave_radiation[j] to
            # the annual-mean insolation of a ZERO-OBLIQUITY planet, Q(phi) = S/pi*cos(phi),
            # constrained so the cos(latitude)-weighted mean is exactly S/4. RMS error
            # 6 W/m2. Hadean obliquity is unknown; zero is the assumption, and it is also
            # the shape the parabola can actually represent — at Earth's 23.44 deg the true
            # curve flattens toward the pole (122 W/m2 there) in a way a parabola cannot
            # follow, giving twice the RMS error. ASSUMPTION — see README.
            ('rad_equator_short', 'ATHAD: TOA short wave insolation at the equator in W/m2', 'double', 298.0),
            ('rad_pole_short', 'ATHAD: TOA short wave insolation at the poles in W/m2', 'double', 0.0),

            ('sigma', 'Stefan-Boltzmann constant W/(m²*K4)', 'double', 5.670280e-8),

            # ==================================================================
            # ATHAD longwave opacity — grey mass absorption coefficients [m2/kg].
            #
            # The layer optical depth is tau = SUM_s kappa_s * q_s * (dp/g) * (p/p_ref),
            # i.e. absorber mass times a pressure-broadening factor. These replace the
            # Bignami/Atwater-Ball emissivity regressions, which are fits to present-day
            # terrestrial columns and saturate to 1 the moment the water path exceeds an
            # Earth-like value.
            #
            # THESE THREE NUMBERS ARE THE BIGGEST SINGLE LEVER ON THE ANSWER and carry
            # roughly a factor-of-two uncertainty. kappa_H2O is the value conventionally
            # used in grey runaway-greenhouse models; CO2 absorbs less per unit mass; the
            # N2/CO/CH4/H2 background is nearly transparent in the thermal infrared. The
            # test of whether they are right is the outgoing longwave flux: a runaway
            # steam atmosphere should sit near the Nakajima / Komabayashi-Ingersoll limit
            # of ~280-310 W/m2, NOT at sigma*T_surf^4.
            ('kappa_H2O', 'ATHAD: grey longwave mass absorption of water vapour in m2/kg', 'double', 0.01),
            ('kappa_CO2', 'ATHAD: grey longwave mass absorption of CO2 in m2/kg', 'double', 0.001),
            ('kappa_bg', 'ATHAD: grey longwave mass absorption of the background gases in m2/kg; LUMPED fallback, used only when ATM_BG_LUMPED=1 (see the per-species values below)', 'double', 1.0e-6),

            # PER-SPECIES BACKGROUND OPACITIES. The single kappa_bg above was 1e-6 — the value
            # of a radiatively inert diatomic — applied to all six background gases at once.
            # Two of them are NOT inert: NH3 and CH4 are 1.4 mole-% each and are among the
            # strongest infrared absorbers in the set, and SO2 is another 1.4 % with strong
            # bands at 7.3 and 8.7 um. Because the six are well mixed and source-free their
            # mass ratios are FIXED, so splitting them collapses to one composition-weighted
            # effective kappa — the radiative transfer does not change shape, only its
            # coefficient. At the configured composition that coefficient is 1.87e-3, i.e.
            # 1867x the lumped value and nearly 2x kappa_CO2, and NH3 (49 %), SO2 (37 %) and
            # CH4 (14 %) carry essentially all of it while N2, H2 and CO contribute ~1 %.
            #
            # THESE ARE ASSUMPTIONS OF THE SAME STANDING AS kappa_H2O AND kappa_CO2 — grey
            # band-averages chosen for the right ORDER and the right ORDERING, not measured.
            # They are set relative to this scheme's own two anchors (H2O 1e-2, CO2 1e-3):
            #   NH3  1e-2  strong 10.5 um band, comparable to water
            #   CH4  3e-3  strong 7.7 um band, between CO2 and water
            #   SO2  2e-3  strong 7.3/8.7 um bands
            #   CO   1e-4  weak dipole, one band at 4.7 um
            #   H2   1e-5  no dipole; collision-induced only, which a grey scheme cannot carry
            #   N2   1e-6  the inert baseline, i.e. what the lumped value assumed for all six
            # A grey scheme cannot represent the window regions these gases fill differently,
            # so the SPLIT is better founded than any individual number in it.
            ('kappa_N2',  'ATHAD: grey longwave mass absorption of N2 in m2/kg',  'double', 1.0e-6),
            ('kappa_CH4', 'ATHAD: grey longwave mass absorption of CH4 in m2/kg', 'double', 3.0e-3),
            ('kappa_NH3', 'ATHAD: grey longwave mass absorption of NH3 in m2/kg', 'double', 1.0e-2),
            ('kappa_H2',  'ATHAD: grey longwave mass absorption of H2 in m2/kg',  'double', 1.0e-5),
            ('kappa_CO',  'ATHAD: grey longwave mass absorption of CO in m2/kg',  'double', 1.0e-4),
            ('kappa_SO2', 'ATHAD: grey longwave mass absorption of SO2 in m2/kg', 'double', 2.0e-3),

            # ATHAD: geothermal / magma-ocean heat flux through the base of the atmosphere
            # in W/m2. A quenching magma ocean radiates far more than the modern Earth's
            # 0.09 W/m2, and at 1500 K this term is plausibly comparable to the absorbed
            # solar — it cannot be omitted. ASSUMPTION.
            ('geothermal_flux', 'ATHAD: heat flux from the molten surface into the atmosphere in W/m2', 'double', 150.0),

            # ATHAD: near-surface Rayleigh (boundary-layer) drag, made CONFIGURABLE by
            # README item 44. Both were constexpr in RHS_Atm_Turb.cpp with comments justifying
            # them by EARTH'S GEOGRAPHY, which is the pattern this project keeps finding:
            #
            #   rayleigh_kf   "baseline 1/day gave ~34 m/s eastward w off W-coast S-America;
            #                  10x cut the surface to ~28 m/s ... so 10x is the settled strength"
            #   drag_n_layers "that ~27 m/s max at 27S/71W is an Andes orographic feature"
            #
            # ATHAD has no South America and no Andes — invariant 1 makes is_land() false
            # everywhere. Worse, drag_n_layers is a COUNT OF CELLS, not a length, so its physical
            # depth is whatever the grid makes it: 5 cells is 236 m on ATOM_Precipitation's grid
            # and 7152 m here, 30.3x deeper. Same defect shape as init_tropopause_layers'
            # round(h / L_atm) — a grid index used as a physical length.
            #
            # DEFAULTS REPRODUCE THE OLD constexpr EXACTLY, so this change alone is
            # bit-identical. They exist so the drag can be SCANNED: the tropical cell decays
            # with no turnover (items 28, 31, 37), drag is the obvious sink, and unlike the
            # radiative questions this one is blocked by neither t_skin (item 43) nor the thermal
            # timescale (item 47) — Psi responds to drag directly and ubud_*/vbud_* can attribute
            # it. Reformulating the depth as a LENGTH in metres is the follow-up; it is kept as a
            # cell count here so that the default stays exactly the measured baseline.
            ('rayleigh_kf', 'ATHAD: near-surface Rayleigh drag rate in 1/s (was a constexpr 1/86400; item 44)', 'double', 1.0/86400.0),
            ('drag_n_layers', 'ATHAD: boundary-layer drag depth in AIR CELLS, not metres (was a constexpr 5.0; item 44 — 5 cells is 7152 m here against 236 m on Earth)', 'double', 5.0),

            ('eps_residuum', 'relative error, end of iterations reached, 1% error  allowed', 'double', 1.0e-4),

#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'none'),
            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_omega_SST'),
#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_omega'),
#            ('turb_model', 'turbulence model: none, k_epsilon, k_omega, k_omega_SST', 'string', 'k_epsilon'),

#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 100),
#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 40),
#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 80),
#            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 300),
            ('inviscid_spinup_iters', 'cumulative iterations to run inviscid (Euler + free-slip mountains) before viscous physics activates; 0 disables', 'int', 0),
            ('inviscid_ramp_iters', 'iterations over which diffusion coefficient ramps from 0 to 1 after the inviscid phase', 'int', 20),

            ('moist_phys_start_iter', 'cumulative iterations before moist physics (SaturationAdjustment, ice scheme, MoistConvection) activates; lets the velocity circulation form on a dry field first; 0 disables (always on)', 'int', 300),

            ('checkpoint_save_iter', 'dump the full 3D prognostic state to output_path/atm_restart_<iter>.bin when total_iter_count reaches this, for a fast debug restart; -1 disables', 'int', 300),
#            ('checkpoint_save_iter', 'dump the full 3D prognostic state to output_path/atm_restart_<iter>.bin when total_iter_count reaches this, for a fast debug restart; -1 disables', 'int', 200),
            ('restart_from_iter', 'load output_path/atm_restart_<iter>.bin and resume from it, skipping the dry spin-up (debug shortcut); -1 disables', 'int', -1),

            # ATHAD: iterations between full 3D restart dumps. Was a constexpr 100 buried in
            # the iteration loop, next to two configurable siblings (checkpoint for the vtk
            # slices, panorama_print for the vts panoramas) — so the one output you cannot
            # switch off was the one that writes 925 MB per dump. 0 or negative disables.
            ('restart_stride', 'ATHAD: iterations between full 3D restart dumps (.bin, ~925 MB each); 0 disables', 'int', 100),
#            ('restart_from_iter', 'load output_path/atm_restart_<iter>.bin and resume from it, skipping the dry spin-up (debug shortcut); -1 disables', 'int', 300),

#            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.001),
#            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.0005),
            # The stretched radial grid is fine at the surface (dr ~ 0.014 vs the old
            # constant 0.025); once the radial finite differences correctly use that
            # spacing, the explicit diffusion CFL limit at the surface tightens to
            # dt ~ dr^2/(2D) ~ 1e-4. dt_visc=5e-4 was ~5x over it and blew up the
            # near-surface cells at iter 174. 1e-4 is CFL-safe (validated: passes 174).
            # ATHAD: im 41 -> 61 shrank the radial step dr from 1/40 to 1/60, and the explicit
            # diffusion CFL limit goes as dr^2 — a factor (40/60)^2 = 0.44 — so this was cut
            # to 4e-5 for the 61-level grid.
            #
            # BACK TO 1e-4 (README items 24, 27), because im = 41 is the default again and
            # the limit scales with the PHYSICAL surface spacing, 1.22 km here against 0.81 km
            # at 61 levels. Not restored on that argument alone: two 200-iteration runs at
            # im = 41, one at each step, both completed and both cleared iteration 174 — where
            # 5e-4 is documented above to have blown up the near-surface cells. Compared at
            # MATCHED PHYSICAL TIME (1e-4 at iteration 80 against 4e-5 at iteration 200, both
            # t = 8.0e-3) the OLR agrees to 1.0 %, 283.55 against 280.75 W/m2, for 2.5x fewer
            # iterations. With im = 41's 1.44x per iteration that is ~3.6x less wall clock to
            # a given physical state, which is what matters against item 18's ~1e4-iteration
            # estimate for geostrophic adjustment.
            #
            # IT IS TIED TO THE GRID. At im = 61 this same value is ~2.2x over the limit; if
            # im goes back up, this must come back down. It was also confirmed with the moist
            # physics running from iteration 0 (item 27, 200 iterations, no NaN), which the
            # dry item-24 scan had not covered.
            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.0001),
            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.000008),
#            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.0005),
#            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.0003),
        ],


        'atmosphere': [

#            ('nm', 'the maximum number of iterations', 'int', 4),
#            ('nm', 'the maximum number of iterations', 'int', 100),
            ('nm', 'the maximum number of iterations', 'int', 400),
            ('checkpoint', "control when to write output files", 'int', 20),

            # Cadence of the TEXT diagnostics (column profile, level summary, OLR), separate
            # from `checkpoint`, which writes ParaView files. 0 = automatic: every 10
            # iterations for a short run (nm <= 100), every 100 for a longer one.
            ('diagnostic_stride', 'ATHAD: iterations between text diagnostics; 0 = auto (10 short / 100 long)', 'int', 0),
            ('panorama_print', "control when to write panorama files", 'int', 100),


            ('coeff_Dalton', "diffusion coefficient in evaporation by Dalton", 'double', 0.7),

#            ('convection_perturbation', 'convective trigger perturbation: 0=fixed, 1=Bechtold (2008) surface-flux-based for shallow/fixed for deep', 'int', 0),
            ('convection_perturbation', 'convective trigger perturbation: 0=fixed, 1=Bechtold (2008) surface-flux-based for shallow/fixed for deep', 'int', 1),

#            ('convection_mode', 'convection type: 0=deep only (precipitating), 1=deep+shallow (non-precipitating if p_diff<p_stat_diff), 2=deep+shallow+midlevel (also non-precipitating for cloud base above p_stat_midlevel/700 hPa)', 'int', 0),
            ('convection_mode', 'convection type: 0=deep only (precipitating), 1=deep+shallow (non-precipitating if p_diff<p_stat_diff), 2=deep+shallow+midlevel (also non-precipitating for cloud base above p_stat_midlevel/700 hPa)', 'int', 1),

#            ('iter_prec', 'precipitation sub-iteration count: min 3 for evaporation (e_d, e_p) to act on non-zero P_conv; check convergence at 4-5', 'int', 4),
            ('iter_prec', 'precipitation sub-iteration count: min 3 for evaporation (e_d, e_p) to act on non-zero P_conv; check convergence at 4-5', 'int', 3),

            ('evap_model', "evaporation formula driving surface humidity update: Dalton, Meyer, or Rohwer", 'string', 'Meyer'),



            # ==================================================================
            # ATHAD vertical grid.
            #
            # L_atm is the AMPLITUDE of the exponential stretch, NOT the shell thickness
            # and NOT a layer spacing. The shell is (exp(zeta) - 1) * L_atm:
            #     Earth: (exp(3.715) - 1) *   400.0 =  16.0 km
            #     ATHAD: (exp(3.000) - 1) * 15719.0 = 300.0 km
            #
            # The shell is set by where the column reaches the radiating level, and the
            # profile has moved under it twice: 300 km on the inherited COSMO shape, 230 km
            # once the proper dry adiabat turned out to compress the atmosphere far more.
            # Finishing the saturation conversion moved it back — with the upper column
            # staying steam it keeps steam's high cp and cools along a shallower adiabat, so
            # 230 km now tops out at 0.29 bar, three times ABOVE the 0.1 bar radiating level,
            # with a top layer of optical depth ~6. An opaque lid emits sigma*T_lid^4 straight
            # out of the domain, which is why the model's OLR is not an output (see t_skin).
            #
            # Deepening it used to be impossible: at 260 km and 300 km the old tridiagonal
            # radiation solve produced NaN across the whole field in its first call, its rows
            # degenerating as eps -> 0. That is fixed (MultiLayerRadiation is now two-stream
            # flux sweeps), and 300 km is what the fix buys:
            #
            #     shell   top p      lid eps   OLR      sigma*T_lid^4
            #     230 km  0.29  bar  1.0000    795 W/m2  787 W/m2   <- OLR == lid, an input
            #     260 km  0.017 bar  0.0610    519       271
            #     300 km  3.8e-4 bar 0.0000    581       271        <- OLR is a real integral
            #
            # At 300 km the top layer is transparent, the isothermal skin is resolved from
            # 256 km up, and the outgoing flux is no longer the boundary temperature read
            # back out. The OLR still varies with the shell depth (519 vs 581), because im is
            # fixed at 61 and a deeper shell is a coarser grid — it is not yet converged.
            #
            # zeta 3.0 with im = 61 keeps the top cell at ~1.7 local scale heights.
            ('L_atm', 'ATHAD: amplitude of the radial stretch in m; shell = (exp(zeta)-1)*L_atm = 300 km', 'double', 15719.0),
            ('zeta', 'ATHAD: radial coordinate-stretching factor (was a hard-coded 3.715)', 'double', 3.0),

            # ATHAD: the radiative-convective boundary of a runaway steam atmosphere sits
            # far higher than Earth's. At 250 km the column is still at 0.68 bar and at
            # 290 km at 0.086 bar, so nearly the whole shell convects. ASSUMPTION — these
            # should be derived from the lapse rate once the radiation is right (Phase 5),
            # not prescribed.
            # ATHAD: the convective column now ends where the adiabat meets the skin
            # temperature, ~207 km. Above that the atmosphere is isothermal and radiative.
            ('tropopause_pole', 'ATHAD: top of the convective column at the poles in m', 'double', 195000.0),
            ('tropopause_equator', 'ATHAD: top of the convective column at the equator in m', 'double', 207000.0),


            # ATHAD: VelocityInitializer prescribes its circulation cells at EARTH latitudes
            # — Hadley 15 deg, Ferrel 45, polar 75, with the trade/westerly nodes between.
            # Those latitudes follow from Earth's thermal Rossby number and are not a
            # property of rotating atmospheres in general:
            #
            #     Ro_T = g*H*(dtheta/theta) / (omega^2 * a^2)
            #
            #             scale height   dtheta/theta   omega      Ro_T     Held-Hou edge
            #     Earth    8.4 km        45/288 = 0.156 7.29e-5    0.0598   18.1 deg
            #     ATHAD   59.3 km        50/1500 = 0.033 3.17e-4   0.0048    5.1 deg
            #
            # Rotation is 4.35x faster (omega^2 18.9x) and the FRACTIONAL equator-pole
            # contrast 4.7x weaker, only partly offset by a 7x deeper atmosphere. A hot
            # surface is not a strongly DIFFERENTIALLY heated one, which is why the higher
            # energy content narrows the circulation instead of widening it — the
            # giant-planet direction rather than the Venus one. The Rhines scale agrees:
            # ~2450 km against Earth's ~3490, so ~8 bands pole-to-pole against ~6.
            #
            # 5.1/15 = 0.34, hence the default 0.33 (Hadley anchor at 5 deg). Four
            # 200-iteration runs, mode 2 balance, moist physics from iteration 0, measured
            # the tropical cell decay against this scale (README item 31):
            #
            #     scale  anchor        Psi @20    Psi @200   decay
            #     1.00   15 deg (Earth) 156 059    138 530   -11.2 %
            #     0.75   11 deg         148 752    134 937    -9.3 %
            #     0.50    7 deg         132 100    122 337    -7.4 %
            #     0.33    5 deg         106 309     98 519    -7.3 %
            #
            # The decay falls monotonically as the cell moves to where the regime wants it
            # and SATURATES between 7 and 5 deg, where Held-Hou puts the edge. The
            # measurement corroborates 0.33 but does not discriminate it from 0.50; it is
            # the Ro_T argument that picks the value, and the two agree.
            #
            # The 50 K contrast is itself prescribed, but the conclusion is not delicate:
            # recovering Earth's Ro_T at this rotation rate would need dT ~ 630 K, and even
            # 200 K still lands near 10 deg.
            #
            # Only the HADLEY anchor is a regime question. Earth's Ferrel cell is thermally
            # indirect and eddy-driven; this flow is axisymmetric to 2 % and neutrally
            # stratified by construction (cosmo_lapse_fraction = 1.0, N^2 ~ 0), so no
            # baroclinic eddies exist to drive one in ANY eon in this model. Scaling those
            # two anchors is meaningless rather than wrong.
            #
            # 1.0 restores Earth's latitudes, which is what every run before README item 31
            # used. The poles are fixed points of the map, so the interpolation still spans
            # the full [0, jm-1].
            ('cell_lat_scale', 'ATHAD: prescribed circulation-cell latitudes as a fraction of Earth\'s; 0.33 = Held-Hou edge for this rotation rate, 1.0 = Earth', 'double', 0.33),

            # ATHAD: cell_lat_scale moves the ANCHORS and nothing else. Every prescribed
            # velocity amplitude in VelocityInitializer stays at the value Earth's bands
            # were tuned to, so the same velocity change is interpolated across a third of
            # the latitude span and EVERY MERIDIONAL GRADIENT IN THE INITIAL STATE IS 1/s
            # TIMES EARTH'S — 3x at the 0.33 default. That is not what "the cell is a third
            # as wide" should mean, and the implied thermal wind goes with it.
            #
            # The scaling laws, taking the latitude scale s = cell_lat_scale and noting
            # that all three modes are a no-op at s = 1:
            #
            #   Continuity.  (1/(a cos)) d(v cos)/dphi + du_r/dz = 0. Shrink the latitude
            #   scale by s and the meridional velocity by s together and dv/dphi is
            #   unchanged, so du_r/dz is unchanged: v scales by s and the RADIAL amplitudes
            #   (ua_00 .. ua_90) do not scale at all. Same vertical motion, narrower cell,
            #   proportionally less meridional flow needed to feed it.
            #
            #   Angular momentum.  A parcel leaving the equator at u = 0 and conserving
            #   angular momentum reaches latitude phi with u = omega*a*sin^2(phi)/cos(phi),
            #   which is proportional to phi^2 for small phi. The jet amplitude ratio
            #   between the compressed and the Earth anchor is therefore s^2 — exact value
            #   at s = 0.33 is 0.110 against s^2 = 0.109, so the small-angle form is good to
            #   1 % here. Earth's observed jets are well below the AMC limit, so this is
            #   applied as a RATIO to the inherited amplitude, not as an absolute jet speed.
            #
            #   0 = no amplitude scaling. What item 32 committed, and what every run before
            #       item 33 used. Preserves the amplitudes and multiplies the shear by 1/s.
            #   1 = v and w by s. Preserves every prescribed meridional gradient exactly:
            #       the initial field is Earth's, geometrically compressed in latitude.
            #   2 = v by s (continuity), w by s^2 (angular momentum). The jet weakens faster
            #       than the overturning, which is what a narrower Hadley cell implies.
            #
            # DEFAULT 0 pending measurement, deliberately: flipping it changes the committed
            # configuration, and this file's own rule is that a default change needs a run
            # behind it. Modes 1 and 2 exist to supply one.
            ('cell_amp_mode', 'ATHAD: scale the prescribed velocity amplitudes with cell_lat_scale; 0 = off (amplitudes unscaled), 1 = v,w by s (continuity), 2 = v by s and w by s^2 (angular momentum)', 'int', 0),

            # ATHAD: how many circulation cells VelocityInitializer lays down per
            # hemisphere. Earth has three -- Hadley, Ferrel, polar -- and the inherited
            # code had that count welded in as ~30 latitude literals plus a hand-written
            # chain of form_diagonals calls. The anchors are now GENERATED:
            #
            #     edges   phi_e(k) = 90*k/n        k = 0..n     (u: ascent/descent branches)
            #     centres phi_c(k) = 90*(k+0.5)/n  k = 0..n-1   (v, w: cell cores)
            #
            # n = 3 returns 0/30/60/90 and 15/45/75, i.e. exactly the inherited anchors.
            #
            # WHY MORE THAN THREE. Earth's three cells follow from Earth's Rhines scale,
            # L_beta = pi*sqrt(2U/beta) with beta = 2*omega*cos(phi)/a. Integrating the
            # band count over the hemisphere, N = a*int_0^(pi/2) sqrt(cos phi) dphi /
            # L_beta(0), gives 2.60 for Earth -- which is the check, since Earth has 3 --
            # and 3.70 here, so ~4. Rotation enters as sqrt(omega) and the wind speed as
            # 1/sqrt(U), so the soft input is U, not the day length: at U = 10 m/s
            # (Earth's) it would be 5.4 cells, at the model's own emergent 21 m/s it is
            # 3.7. Treat "4" as "4 to 5, uncertain in the wind".
            #
            # WHAT THE EXTRA CELL IS. The amplitudes do not generalise on their own,
            # because the inherited polar cell carries the FERREL's sense (v_trop +0.5
            # against Hadley's -3.0), so the direct/indirect alternation that would let
            # any n be tiled is not what the table actually contains. The rule chosen is
            # therefore explicit: cell 0 keeps the Hadley template, cell n-1 keeps the
            # polar template, and cells 1..n-2 are Ferrel copies -- the extra band is
            # inserted in mid-latitudes, where the Rhines argument says the deformation
            # scale shrinks, rather than at the pole. At n = 3 that is Hadley/Ferrel/polar
            # unchanged.
            #
            # DEFAULT 5, changed from Earth's 3 after items 36-37. Three independent
            # arguments, none of them a converged measurement:
            #
            #   1. Band width. With the Hadley edge at 30*s = 9.9 deg, tiling the rest of
            #      the hemisphere gives extratropical cells of 40 deg at n = 3, 27 at
            #      n = 4 and 20 at n = 5, against the Rhines scale's ~22. Only n = 5
            #      matches, and it is the same layout the hemisphere-integrated band count
            #      N = a*int sqrt(cos phi) dphi / L_beta(0) points at once the narrow
            #      direct cell is counted separately (~1 + 4 rather than a uniform ~3.7).
            #   2. Geometry. At n = 3 and s = 0.33 the two extratropical cells are 40 deg
            #      wide -- wider than Earth's 30 -- so compressing the tropics visibly
            #      stretches everything else. n = 5 divides the extratropics evenly.
            #   3. Core position, and this is the only one that is not a rate. At 400
            #      iterations n = 5's Psi maxima still sit exactly on its prescribed cores
            #      at 20/40/60 deg, while n = 3's 30 deg core is no longer a peak at all,
            #      the maximum having migrated equatorward past 20. The model keeps the
            #      layout n = 5 gives it and rearranges the one n = 3 gives it.
            #
            # WHAT THE DEFAULT DOES NOT REST ON. Item 36 also reported that n = 5 deepens
            # its interior boundaries; item 37 withdrew that, because by 400 iterations
            # every boundary in both layouts erodes to the same ~3:1 contrast. And the
            # tropical cell does not care about n at all -- Psi_max decays -16.7 % at n = 3
            # against -16.2 % at n = 5 over 20->400, converging to within 0.35 %.
            #
            # CAVEAT, the same one item 31 raised about the Ferrel and polar anchors: this
            # flow is axisymmetric to 2 % and neutrally stratified by construction, so there
            # are no baroclinic eddies to maintain any indirect cell. n sets the structure
            # the model is handed, not a structure it can generate. Set 3 to recover Earth's
            # layout, which is what everything before item 37 was measured on.
            ('n_cells_hemisphere', 'ATHAD: prescribed circulation cells per hemisphere; 5 = this rotation rate (Rhines ~22 deg bands), 3 = Earth (Hadley/Ferrel/polar)', 'int', 5),


            # ATHAD albedo. These replace the inherited albedo_pole/albedo_equator, which
            # were INERT — MultiLayerRadiation built its own albedo from bare literals and
            # never read them, so the pole/equator pair was configuration theatre.
            #
            # albedo_surface: a quenching silicate melt is dark. Measured basaltic-melt
            # albedos are 0.05-0.10. There is no land, no snow and no sea ice, so this is a
            # single global clear-sky value, not Earth's latitude parabola.
            #
            # albedo_cloud: the SW albedo of an optically thick cloud top, composited over
            # the surface value by refl = tau/(tau+2) on the condensate path. IT IS
            # CURRENTLY THE WHOLE ANSWER. The reflectivity saturates the moment any
            # condensate is present, and ATHAD's cloud deck is enormously thick, so the
            # model's mean planetary albedo IS this number to four decimal places. It is
            # therefore the second-biggest lever after the opacities, and it is an
            # ASSUMPTION: 0.50 is a thick terrestrial water cloud. A deep, cold, slowly
            # sedimenting Hadean deck could plausibly be brighter.
            ('albedo_surface', 'ATHAD: clear-sky albedo of the molten silicate surface', 'double', 0.08),
            ('albedo_cloud', 'ATHAD: shortwave albedo of an optically thick cloud top', 'double', 0.50),

            ('epsilon_equator', 'emissivity and absorptivity caused by other gases than water vapour/(by Häckel)', 'double', 0.48),
            ('epsilon_pole', 'emissivity and absorptivity caused by other gases than water vapour at the poles', 'double', 0.45),
            ('epsilon_tropopause', 'emissivity and absorptivity caused by other gases than water vapour in the tropopause', 'double', 0.001),

            ('re', 'Reynolds number for laminar flows: ratio viscous to inertia forces, Re = u * L/nue', 'double', 1000.0),
            ('sc_WaterVapour', 'Schmidt number of water vapour, Sc = nue/D', 'double', 0.61),
            ('sc_CO2', 'Schmidt number of CO2', 'double', 0.96),
            ('pr', 'Prandtl number of air for laminar flows', 'double', 0.7179),
            ('pr_turb', 'turbulent Prandtl number for temperature transport in turbulent flows', 'double', 0.9),
            # ATHAD: the boundary layer scales with the scale height, which is 59 km here
            # against Earth's 8.4 km — a ~7x ratio applied to Earth's 1500 m. At the chosen
            # grid this puts ~12 cells inside the ABL. ASSUMPTION.
            ('abl_height', 'ATHAD: physical depth of the atmospheric boundary layer in m', 'double', 10000.0),
            # ==================================================================
            # ATHAD atmospheric composition — MOLE fractions of the Hadean mixture.
            # The residual 7% beyond H2O/CO2/N2 is split evenly across five trace gases.
            # These are the INPUT; MixtureAtm.h derives the mass fractions and the
            # mixture gas constant from them and checks they sum to 1.
            # ==================================================================
            ('x_H2O', 'ATHAD: mole fraction of H2O', 'double', 0.800),
            ('x_CO2', 'ATHAD: mole fraction of CO2', 'double', 0.100),
            ('x_N2',  'ATHAD: mole fraction of N2',  'double', 0.030),
            ('x_CH4', 'ATHAD: mole fraction of CH4', 'double', 0.014),
            ('x_NH3', 'ATHAD: mole fraction of NH3', 'double', 0.014),
            ('x_H2',  'ATHAD: mole fraction of H2',  'double', 0.014),
            ('x_CO',  'ATHAD: mole fraction of CO',  'double', 0.014),
            ('x_SO2', 'ATHAD: mole fraction of SO2', 'double', 0.014),

            # ep = R_background / R_H2O = 317.3/461.5. NOTE: the dilute approximation this
            # constant serves, q_sat = ep*E/(p-(1-ep)*E), is INVALID here because H2O is the
            # bulk gas, not a trace. Phase 4 replaces it with the exact mass-fraction form;
            # ep remains only where a genuine gas-constant ratio is wanted.
            ('ep', 'ATHAD: ratio of the background-mixture to water-vapour gas constants', 'double', 0.6875),
            ('hp', 'water vapour pressure at T = 0°C: E = 6.1 hPa', 'double', 6.1078),

            # ATHAD: "Air" now means the NON-CONDENSABLE BACKGROUND (everything but H2O and
            # CO2): x_bg = 0.100, M_bg = 26.207 g/mol -> R_bg = 317.3 J/(kg K). It is not air.
            ('R_Air', 'ATHAD: specific gas constant of the non-condensable background in J/(kg*K)', 'double', 317.3),
            ('R_WaterVapour', 'specific gas constant of water vapour in J/(kg*K)', 'double', 461.5),
            ('r_water_vapour', 'density of saturated water vapour in kg/m³ at 10°C', 'double', 0.0094),
            ('R_co2', 'specific gas constant of CO2 in J/(kg*K)', 'double', 188.9),
            ('lv', 'specific latent evaporation heat(condensation heat) in J/kg', 'double', 2.52e6),
            ('ls', 'specific latent vaporisation heat(sublimation heat) in J/kg', 'double', 2.83e6),

            # ATHAD: cp of the MIXTURE at Hadean temperatures, ~2x Earth's 1005. Mass-weighted
            # from H2O ~2400, CO2 ~1280, background ~1300 J/(kg K) at 1000-1500 K. This constant
            # is the fallback; MixtureAtm::cp_of() gives the local, temperature-dependent value.
            # cv_l = cp_l - R_mix = 2040 - 387.9.
            ('cp_l', 'ATHAD: specific heat capacity of the mixture at constant pressure in J/(kg K)', 'double', 2040.0),
            ('cv_l', 'ATHAD: specific heat capacity of the mixture at constant volume in J/(kg K)', 'double', 1652.1),
            ('lamda', 'heat transfer coefficient of air in W/(m K)', 'double', 0.0262),
            ('r_co2', 'density of CO2 in kg/m³ at 25°C', 'double', 0.0019767),
            ('gam', 'constant slope of temperature    gam = 6.5 K/1000 m', 'double', 0.0065),

            ('u_0', 'annual mean of surface wind velocity in m/s, 8 m/s compare to 28.8 km/h', 'double', 8.0),
            # ATHAD: an upper PHYSICAL bound on the prognostic temperature, replacing the
            # hard-coded 333.15 K (60 °C) literal in SaturationAdjustment. That literal was
            # justified as "well above any physical surface temperature" — true on Earth,
            # but a factor of 4.5 BELOW ATHAD's 1500 K surface, and it was written back into
            # the prognostic field, collapsing 250 bar to 30 bar in one iteration.
            # 2000 K leaves headroom over the surface and stays inside the cp Shomate fits.
            # ATHAD: ceiling on the condensate mass fraction, applied in two places that
            # were both bare literals — SaturationAdjustment::clampAndFade caps cloud and
            # ice SEPARATELY at this value, and the moist-physics block in RunTimeSlice caps
            # their SUM at it before densities() reads them. Kept as one number so they
            # cannot drift apart; the sum-cap is the binding one.
            #
            # The inherited comment called 0.05 "~50x the largest physical cloud/ice mixing
            # ratio, so it never clips a real cloud — it only stops a runaway". That is an
            # EARTH statement: terrestrial cloud water is ~1 g/kg. ATHAD's atmosphere is
            # 67 % water by mass and its polar column pegs this cap exactly (49.999996 g/kg
            # measured at 80N, 185.6 km), so here it is not a runaway backstop at all — it
            # is setting the condensate, and with it the optical depth and the albedo.
            # Raise it and re-measure before trusting a cloud field that sits on it.
            ('cloud_cap', 'ATHAD: ceiling on the condensate mass fraction (cloud, ice, and their sum) in kg/kg', 'double', 0.05),

            ('t_max_phys', 'ATHAD: upper physical bound on the prognostic temperature in K', 'double', 2000.0),

            ('t_00', 'temperature in K compare to -37°C', 'double', 236.15),
            ('t_000', 'temperature in K compare to -20°C', 'double', 253.15),
            ('s_0', 'entropy at 0°C, cp_l * t_0 in J/kg', 'double', 274515.75),
            # ATHAD: the water-vapour scale is the Hadean mass fraction q_H2O = 0.6724, not
            # Earth's 0.035 trace. c_0 is a normalisation in the RHS energy/moisture
            # coefficients (RHS_Atm_Turb.cpp: coeff_energy, coeff_MC_q, coeff_L).
            ('c_0', 'ATHAD: reference water vapour mass fraction in kg/kg', 'double', 0.6724),

            # ATHAD: the CO2 field is a MASS FRACTION, not ppm. At 20.5% by mass, ppm is
            # meaningless. The Hadean has no biosphere, no vegetation and no carbonate ocean
            # sink, so CO2 is simply well mixed — the Earth surface-source/tropopause-sink
            # parabola and the vegetation/ocean/land ppm budgets have no subject here.
            ('co2_0', 'ATHAD: reference CO2 mass fraction in kg/kg', 'double', 0.2053),
            ('co2_scale', 'multiplier applied to the whole CO2 field for sensitivity experiments (1.0 = field as built; 2.0 = doubled CO2)', 'double', 1.0),

            # ATHAD: no land, so no land/ocean humidity split — a single surface relative
            # humidity applies everywhere. Kept as one knob rather than two identical ones.
            ('c_ocean', 'ATHAD: surface water vapour as a fraction of the saturation value, in %', 'double', 70),

        ],
    }




    XML_READ_FUNCS = {                                                  # dictionary (XML_READ_FUNCS) with keys ('"string", etc.) and one list element ("FillStringWithElement")
        "string": "FillStringWithElement",
        "double": "FillDoubleWithElement",
        "int": "FillIntWithElement",
        "bool": "FillBoolWithElement"
    }



    # functions begin

    def write_cpp_defaults(filename, classname, sections):

        with open(filename, 'w') as f:

            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")
            f.write("void %s::SetDefaultConfig() {\n" % classname)

            for section in sections:
                f.write('\n  // %s section\n' % section)

                for slug, desc, ctype, default in PARAMS[section]:
                    rhs = default

                    if ctype == 'string':
                        rhs = '"%s"' % default
                    elif ctype == 'bool':

                        if default:
                            rhs = 'true'
                        else:
                            rhs = 'false'

                    f.write('  %s = %s;\n' %(slug, rhs))

            f.write("}")



    def write_cpp_load_config(filename, classname, sections):

        with open(filename, 'w') as f:

            f.write("// config files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")

            for section in sections:
                f.write('\n  // %s section\n' % section)
                element_var_name = 'elem_%s' % section
                f.write('\n  if(%s) {\n' %(element_var_name))

                for slug, desc, ctype, default in PARAMS [section]:
                    func_name = XML_READ_FUNCS [ctype]
                    f.write('    Config::%s(%s, "%s", %s);\n' %(func_name, element_var_name, slug, slug))

                f.write("  }\n")



    def write_cpp_headers(filename, sections, is_extern = False):

        with open(filename, 'w') as f:

            f.write("// header files\n")
            f.write("// THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("// ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write("\n")

            if is_extern:
                f.write("#include<string>\n\n")
                f.write("using namespace std;\n")

                #if 'atmosphere' in filename:
                #    f.write("namespace AtmParameters{\n")
                #else:
                #    f.write("namespace HydParameters{\n")

            for section in sections:
                f.write('\n// %s section\n' % section)

                for slug, desc, ctype, default in PARAMS [section]:
                    if is_extern:
                        f.write('   extern %s %s;\n' %(ctype, slug))
                    else:
                        f.write('%s %s;\n' %(ctype, slug))
           
            if is_extern:
                f.write("}\n")



    def write_pxi(input_filename, output_filename, substitutions):

        data = open(input_filename, 'r').read()
        indent = '    '

        for key, classname, sections in substitutions:
            rep = ''

            for section in sections:
                rep += '%s# %s section\n' %(indent, section)

                for slug, desc, ctype, default in PARAMS[section]:
                    rep += '%sproperty %s:\n' %(indent, slug)
                    rep += '%s    def __get__(%s self):\n' %(indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        return self._thisptr.%s\n' %(indent, slug)
                    rep += '%s\n' % indent
                    rep += '%s    def __set__(%s self, value):\n' %(indent, classname)
                    rep += '%s        self._check_alive()\n' % indent
                    rep += '%s        self._thisptr.%s = <%s> value\n' %(indent, slug, ctype)
                    rep += '%s\n' % indent

            data = data.replace('{{ %s }}' % key, rep)

        with open(output_filename, 'w') as f:

            f.write("""# pxi files\n""")
            f.write("# THIS FILE IS AUTOMATICALLY GENERATED BY param.py\n")
            f.write("# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME\n")
            f.write(data)



    def write_pxd(filename, model, sections):

        with open(filename, 'w') as f:
            # Sadly, Cython docs are incorrect on usage of 'include', so we must include a whole lot of boilerplate

            f.write("""# pxd files\n""")
            f.write("""# THIS FILE IS AUTOMATICALLY GENERATED BY param.py
# ANY CHANGES WILL BE OVERWRITTEN AT COMPILE TIME
from libcpp.vector cimport vector
cdef extern from "c%sModel.h":
    cppclass c%sModel:
        c%sModel() except +  # NB! std::bad_alloc will be converted to MemoryError
        void LoadConfig(const char *filename)
        void Run()
        void RunTimeSlice(int time_slice)
        vector[float] get_layer_heights()
""" %(model, model, model))

            for section in sections:
                f.write('        # %s section\n' % section)

                for slug, desc, ctype, default in PARAMS [section]:
                    f.write('        %s %s\n' %(ctype, slug))



    def write_config_xml(filename, sections):

        with open(filename, 'w') as f:
            f.write("""<!-- THIS FILE IS GENERATED AUTOMATICALLY BY param.py. DO NOT EDIT. -->""")
            f.write('<atom>')

            for section in sections:
                f.write('    <%s>\n' % section)

                for slug, desc, ctype, default in PARAMS [section]:
                    if ctype == 'bool':
                        default = str(default).lower()                  # Python uses True/False, C++, uses true/false
                    f.write('        <%s>%s</%s>  <!-- %s(%s) -->\n' %(slug, default, slug, desc, ctype))

                f.write('    </%s>\n' % section)
 
            f.write('</atom>')

    # functions end



    atmosphere_sections = ['common', 'atmosphere']

    for filename, classname, sections in [
        ('atmosphere/cAtmosphereDefaults.cpp.inc', 'cAtmosphereModel', atmosphere_sections)
    ]:
        write_cpp_defaults(filename, classname, sections)


    for filename, classname, sections in [
        ('atmosphere/AtmosphereLoadConfig.cpp.inc', 'cAtmosphereModel', atmosphere_sections)
    ]:
        write_cpp_load_config(filename, classname, sections)


    for filename, sections in [
        ('atmosphere/AtmosphereParams.h.inc', atmosphere_sections)
    ]:
        write_cpp_headers(filename, sections)


    write_pxi('python/pyathad.pyx.template', 'python/pyathad.pyx', [
        ('atmosphere_params', 'Atmosphere', atmosphere_sections)])


    for filename, model, sections in [
        ('python/atmosphere_pxd.pxi', 'Atmosphere', atmosphere_sections)
    ]:
        write_pxd(filename, model, sections)


    for  filename, sections in [
        ('python/config_athad.xml', atmosphere_sections)
    ]:
        write_config_xml(filename, sections)


    for  filename, sections in [
        ('cli/config_athad.xml', atmosphere_sections)
    ]:
        write_config_xml(filename, sections)



if __name__ == '__main__':
    main()
