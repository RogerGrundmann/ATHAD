# Given a parameter definition, generates necessary C++, Python and XML bindings
# coding=utf-8



def main():

    # read the input definition
    # name, description, datatype, default in the tuples


    PARAMS = {                                                          # dictionary{} (PARAMS) with keys ('common', etc.):[ and their tuples ('output_path', etc.)]
        'common': [
            ('output_path', 'directory where model outputs should be placed(must end in /)', 'string', 'output-Hadean/'),

            ('bathymetry_path', 'directory where the topografic grids are located', 'string', '../data/topo_grids'),
#            ('bathymetry_path', 'directory where the topografic grids are located', 'string', '../data/simon_topo'),

            ('BathymetrySuffix', 'suffix of the timesteps Ma in million years', 'string', 'Ma_smooth.xyz'),
#            ('BathymetrySuffix', 'suffix of the timesteps Ma in million years', 'string', 'Ma_Simon.xyz'),

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
            ('omega', 'rotation rate of the earth in rad/s', 'double', 7.292e-5),

            ('g', 'gravitational acceleration of the earth in m/s²', 'double', 9.8066),








            #parameters for data reconstruction

            ('time_start', 'start time', 'int', 0),
#            ('time_end', 'end time', 'int', 10),
            ('time_end', 'end time', 'int', 0),
            ('time_step', 'step size between timeslices', 'int', 10),

            ('velocity_v_file', '' ,'string','../data/v_surface.txt'),
            ('velocity_w_file', '' ,'string','../data/w_surface.txt'),

            ('temperature_file', '', 'string', '../data/SurfaceTemperature_NASA.xyz'),

            ('precipitation_file', '', 'string', '../data/SurfacePrecipitation_NASA.xyz'),

            ('salinity_file', '', 'string', '../data/SurfaceSalinity_NASA.xyz'),

            ('temperature_global_file', '', 'string', '../data/scotese_etal_2021_global_temp_1my.txt'),
            ('temperature_equat_file', '', 'string', '../data/scotese_etal_2021_equat_temp_1my.txt'),
            ('temperature_pole_file', '', 'string', '../data/scotese_etal_2021_polar_temp_1my.txt'),

            ('reconstruction_script_path', '', 'string', '../reconstruction/reconstruct_atom_data.py'),

            ('use_earthbyte_reconstruction', 'control whether use earthbyte method to recontruct grids', 'bool', False),
#            ('use_earthbyte_reconstruction', 'control whether use earthbyte method to recontruct grids', 'bool', True),

            ('use_NASA_velocity', 'if use NASA velocity to initialise velocity', 'bool', False),
#            ('use_NASA_velocity', 'if use NASA velocity to initialise velocity', 'bool', True),

            ('use_NASA_temperature', 'if use NASA temperature to initialise surface temperature', 'bool', True),

            ('use_NASA_salinity', 'if use NASA sea-surface salinity to initialise surface salinity: Ma=0 uses the 2D field, Ma>0 its zonal-mean latitude profile; false = invert the Gill density equation', 'bool', True),

#            ('Ma_switch', 'switch initial temperatur from NASA to parabolic approach', 'int', 50),
            ('Ma_switch', 'switch initial temperatur from NASA to parabolic approach', 'int', 100),

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
            # near-surface lapse rate beta*g/(R*T0). The inherited beta = 42 K is an EARTH
            # constant: it gives 4.99 K/km at R=286.9, T0=288 — about 0.511 of the dry adiabat
            # (beta_adiabatic = R*T0/cp = 82.2), i.e. tuned to the observed moist lapse.
            #
            # Carried over unchanged to ATHAD it would give 0.71 K/km — an essentially
            # ISOTHERMAL 300 km column, because beta*g/(R*T0) falls with both the larger R and
            # the far larger T0. So beta is derived rather than fixed:
            #
            #     beta = cosmo_lapse_fraction * R_mix * T_surf / cp
            #
            # which reproduces Earth's beta at Earth's numbers and gives ~146 K and 2.46 K/km
            # here. A steam atmosphere releasing this much latent heat should indeed sit well
            # below its dry adiabat (4.81 K/km), so a sub-adiabatic fraction is right — but the
            # VALUE 0.511 is inherited Earth tuning, not a Hadean result.
            ('cosmo_lapse_fraction', 'ATHAD: near-surface lapse rate as a fraction of the dry adiabat; sets the COSMO beta', 'double', 0.5108),
            # ATHAD reference density of the MIXTURE at the surface, not of dry air:
            # rho = p/(R_mix*T) = 25e6 Pa / (387.9 * 1500 K) = 42.97 kg/m³ (Earth: 1.2041).
            # This is what sets the surface pressure, via p = 1e-2*(r_air*R_Air*T).
            ('r_air', 'ATHAD: reference density of the atmospheric mixture at the surface in kg/m³', 'double', 42.97),
            ('r_0_water', 'reference density of fresh water in kg/m3', 'double', 997.0),
            ('t_equat_modern', 'mean temperature of the modern earth in °C', 'double', 15.4),
            ('t_pole_modern', 'pole temperature of the modern earth in °C', 'double', - 15.4),

            ('t_paleo_max', 'maximum add of mean temperature in °C during paleo times', 'double', 10.0),
 
            # ATHAD insolation: the faint young Sun. At 4.4 Ga the solar constant was about
            # 0.71 of today's 1361 W/m2, i.e. 966 W/m2. Spread over a rotating sphere the
            # global mean absorbed-at-TOA figure is S/4 = 242 W/m2 before albedo; the
            # equator/pole split below carries the same latitudinal contrast the Earth
            # values did, scaled by 0.71. ASSUMPTION — see README.
            ('rad_equator_short', 'ATHAD: short wave radiation at the equator in W/m2', 'double', 116.0),
            ('rad_pole_short', 'ATHAD: short wave radiation at the poles in W/m2', 'double', 71.0),

            # These two are the LONGWAVE boundary values of the inherited scheme. They are
            # present-day Earth fluxes and have no Hadean meaning; mode 2 computes the
            # longwave from the optical depth instead, so they survive only where the old
            # code still reads them.
            ('rad_equator', 'inherited Earth longwave boundary value in W/m2 (unused in radiation_mode 2)', 'double', 398.2),
            ('rad_pole', 'inherited Earth longwave boundary value in W/m2 (unused in radiation_mode 2)', 'double', 360.0),

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
            ('kappa_bg', 'ATHAD: grey longwave mass absorption of the background gases in m2/kg', 'double', 1.0e-6),

            # ATHAD: geothermal / magma-ocean heat flux through the base of the atmosphere
            # in W/m2. A quenching magma ocean radiates far more than the modern Earth's
            # 0.09 W/m2, and at 1500 K this term is plausibly comparable to the absorbed
            # solar — it cannot be omitted. ASSUMPTION.
            ('geothermal_flux', 'ATHAD: heat flux from the molten surface into the atmosphere in W/m2', 'double', 150.0),

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
#            ('restart_from_iter', 'load output_path/atm_restart_<iter>.bin and resume from it, skipping the dry spin-up (debug shortcut); -1 disables', 'int', 300),

#            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.001),
#            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.0005),
            # The stretched radial grid is fine at the surface (dr ~ 0.014 vs the old
            # constant 0.025); once the radial finite differences correctly use that
            # spacing, the explicit diffusion CFL limit at the surface tightens to
            # dt ~ dr^2/(2D) ~ 1e-4. dt_visc=5e-4 was ~5x over it and blew up the
            # near-surface cells at iter 174. 1e-4 is CFL-safe (validated: passes 174).
            # ATHAD: im 41 -> 61 shrinks the radial step dr from 1/40 to 1/60, and the
            # explicit diffusion CFL limit goes as dr^2 — a factor (40/60)^2 = 0.44. The
            # inherited 1e-4 was validated at im=41 and would be ~2.2x over the limit here,
            # so it is scaled to 4e-5. Verify against a long run before trusting it.
            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.00004),
            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.000008),
#            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.0005),
#            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.0003),
        ],


        'atmosphere': [

#            ('nm', 'the maximum number of iterations', 'int', 4),
#            ('nm', 'the maximum number of iterations', 'int', 100),
            ('nm', 'the maximum number of iterations', 'int', 400),
            ('checkpoint', "control when to write output files", 'int', 20),
            ('panorama_print', "control when to write panorama files", 'int', 100),


            ('coeff_Dalton', "diffusion coefficient in evaporation by Dalton", 'double', 0.7),

#            ('convection_perturbation', 'convective trigger perturbation: 0=fixed, 1=Bechtold (2008) surface-flux-based for shallow/fixed for deep', 'int', 0),
            ('convection_perturbation', 'convective trigger perturbation: 0=fixed, 1=Bechtold (2008) surface-flux-based for shallow/fixed for deep', 'int', 1),

#            ('convection_mode', 'convection type: 0=deep only (precipitating), 1=deep+shallow (non-precipitating if p_diff<p_stat_diff), 2=deep+shallow+midlevel (also non-precipitating for cloud base above p_stat_midlevel/700 hPa)', 'int', 0),
            ('convection_mode', 'convection type: 0=deep only (precipitating), 1=deep+shallow (non-precipitating if p_diff<p_stat_diff), 2=deep+shallow+midlevel (also non-precipitating for cloud base above p_stat_midlevel/700 hPa)', 'int', 1),

#            ('iter_prec', 'precipitation sub-iteration count: min 3 for evaporation (e_d, e_p) to act on non-zero P_conv; check convergence at 4-5', 'int', 4),
            ('iter_prec', 'precipitation sub-iteration count: min 3 for evaporation (e_d, e_p) to act on non-zero P_conv; check convergence at 4-5', 'int', 3),

            ('evap_model', "evaporation formula driving surface humidity update: Dalton, Meyer, or Rohwer", 'string', 'Meyer'),


            ('Ma_max', 'parabolic temperature distribution 300 Ma(from Ruddiman)', 'int', 300),
            ('Ma_max_half', 'half of time scale', 'int', 150),

            # ==================================================================
            # ATHAD vertical grid.
            #
            # L_atm is the AMPLITUDE of the exponential stretch, NOT the shell thickness
            # and NOT a layer spacing. The shell is (exp(zeta) - 1) * L_atm:
            #     Earth: (exp(3.715) - 1) *   400.0 =  16.0 km
            #     ATHAD: (exp(3.000) - 1) * 15718.7 = 300.0 km
            #
            # 300 km is set by where the column reaches the radiating level. With R_mix =
            # 387.9 and T_surf = 1500 K the scale height is 59 km at the surface, and the
            # COSMO profile terminates (T -> 0) at 305 km. Measured top pressures:
            #     150 km -> 13.07 bar     250 km -> 0.679 bar
            #     200 km ->  3.58 bar     300 km -> 0.033 bar   <- below the 0.1 bar target
            #
            # zeta was reduced 3.715 -> 3.0 and im raised 41 -> 61 for resolution ALOFT.
            # The stretch concentrates levels near the surface, where the scale height is
            # largest, and coarsens them aloft where it is smallest — backwards for
            # pressure resolution. Top-cell thickness in local scale heights:
            #     im=41, zeta=3.715 -> 2.92 H   (a cell spanning ~3 scale heights)
            #     im=61, zeta=3.0   -> 1.65 H   <- chosen
            #     im=81, zeta=2.5   -> 1.08 H   (better, at 2x the im=41 cost)
            # Surface spacing at the chosen setting is 806 m (Earth: 39 m), which is fine
            # against a 59 km scale height.
            ('L_atm', 'ATHAD: amplitude of the radial stretch in m; shell = (exp(zeta)-1)*L_atm = 300 km', 'double', 15718.7),
            ('zeta', 'ATHAD: radial coordinate-stretching factor (was a hard-coded 3.715)', 'double', 3.0),

            # ATHAD: the radiative-convective boundary of a runaway steam atmosphere sits
            # far higher than Earth's. At 250 km the column is still at 0.68 bar and at
            # 290 km at 0.086 bar, so nearly the whole shell convects. ASSUMPTION — these
            # should be derived from the lapse rate once the radiation is right (Phase 5),
            # not prescribed.
            ('tropopause_pole', 'ATHAD: extension of the troposphere at the poles in m', 'double', 250000.0),
            ('tropopause_equator', 'ATHAD: extension of the troposphere at the equator in m', 'double', 280000.0),


            # ATHAD: Earth's 0.294/0.1 split encodes polar ice and open ocean, neither of
            # which exists here. A runaway steam atmosphere is expected to carry a thick
            # global cloud deck, so a single high, latitude-independent value is more
            # defensible than a contrast built from surfaces that are not present.
            # ASSUMPTION, and a strong lever on the absorbed solar.
            ('albedo_pole', 'ATHAD: cloud-deck albedo at the poles', 'double', 0.4),
            ('albedo_equator', 'ATHAD: cloud-deck albedo at the equator', 'double', 0.4),

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

            ('sst_coupling_alpha', 'outer-loop (Picard) hydrosphere->atmosphere SST coupling strength: blend fraction of the hydrosphere surface SST (read from <stem>_Transfer_Hyd_SST_<iter>.vwtp) into the atmospheric ocean surface temperature t.x[0] at init, t.x[0] <- (1-alpha)*t.x[0] + alpha*SST_hyd, before the t_eq snapshot so it propagates into the Held-Suarez target. 0.0 = OFF (no read; identical to the one-way chain and to round 0 of a Picard loop, which has no SST file yet). Ocean-only, finite-checked, SST-clamped to [-1.8,40] C, ocean-mean-anchored so total energy does not drift. Under-relax across rounds (e.g. 0.3-0.5)', 'double', 0.0),
            ('hyd_sst_iter', 'which hydrosphere SST snapshot to read for sst_coupling_alpha: reads <stem>_Transfer_Hyd_SST_<iter>.vwtp for this iteration; -1 = use the latest (highest-iter) snapshot present in the output dir', 'int', -1),
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
