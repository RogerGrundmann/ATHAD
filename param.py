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

            ('p_0', 'pressure at sea level in hPa', 'double', 1013.25),
            ('t_0', 'temperature in K compare to 0°C', 'double', 273.15),
            ('r_air', 'density of dry air in kg/m³ at 20°C', 'double', 1.2041),
            ('r_0_water', 'reference density of fresh water in kg/m3', 'double', 997.0),
            ('t_equat_modern', 'mean temperature of the modern earth in °C', 'double', 15.4),
            ('t_pole_modern', 'pole temperature of the modern earth in °C', 'double', - 15.4),

            ('t_paleo_max', 'maximum add of mean temperature in °C during paleo times', 'double', 10.0),
 
            ('rad_equator', 'long wave radiation in W/m2, t_equator = 1.0976 compares to 28.0°C = 299.81 K', 'double', 398.2),
            ('rad_pole', 'long wave radiation in W/m2, t_pole = 0.9436 compares to -15.4°C = 257.75 K', 'double', 360.0),
#            ('rad_equator', 'long wave radiation in W/m2, t_equator = 1.0976 compares to 28.0°C = 299.81 K', 'double', 239.9),
#            ('rad_pole', 'long wave radiation in W/m2, t_pole = 0.9436 compares to -15.4°C = 257.75 K', 'double', 250.26),

            ('rad_equator_short', 'short wave radiation in W/m2', 'double', 163.3),
            ('rad_pole_short', 'short wave radiation in W/m2', 'double', 100.0),

            ('sigma', 'Stefan-Boltzmann constant W/(m²*K4)', 'double', 5.670280e-8),

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
            ('dt_visc', 'non-dimensional time step used in the viscous (production) phase', 'double', 0.0001),
            ('dt_inviscid', 'non-dimensional time step used during the inviscid spin-up phase (smaller to absorb the missing diffusive damping)', 'double', 0.00002),
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

            ('L_atm', 'extension of the atmosphere shell in m, total height is 16000m*40 steps', 'double', 400.0),

            ('tropopause_pole', 'extension of the troposphere at the poles in m', 'double', 8000.0),
            ('tropopause_equator', 'extension of the troposphere at the equator in m', 'double', 15000.0),


            ('albedo_pole', 'albedo around the poles', 'double', 0.294),
            ('albedo_equator', 'albedo around the equator', 'double', 0.1),

            ('epsilon_equator', 'emissivity and absorptivity caused by other gases than water vapour/(by Häckel)', 'double', 0.48),
            ('epsilon_pole', 'emissivity and absorptivity caused by other gases than water vapour at the poles', 'double', 0.45),
            ('epsilon_tropopause', 'emissivity and absorptivity caused by other gases than water vapour in the tropopause', 'double', 0.001),

            ('re', 'Reynolds number for laminar flows: ratio viscous to inertia forces, Re = u * L/nue', 'double', 1000.0),
            ('sc_WaterVapour', 'Schmidt number of water vapour, Sc = nue/D', 'double', 0.61),
            ('sc_CO2', 'Schmidt number of CO2', 'double', 0.96),
            ('pr', 'Prandtl number of air for laminar flows', 'double', 0.7179),
            ('pr_turb', 'turbulent Prandtl number for temperature transport in turbulent flows', 'double', 0.9),
            ('abl_height', 'physical depth of the atmospheric boundary layer in m, sets the top of the surface-driven turbulent TKE seeding profile (decoupled from the L_atm grid length scale)', 'double', 1500.0),
            ('ep', 'ratio of the gas constants of dry air to water vapour [kg_air/kg_vapour]', 'double', 0.62198),
            ('hp', 'water vapour pressure at T = 0°C: E = 6.1 hPa', 'double', 6.1078),
            ('R_Air', 'specific gas constant of air in J/(kg*K)', 'double', 286.9),
            ('R_WaterVapour', 'specific gas constant of water vapour in J/(kg*K)', 'double', 461.4),
            ('r_water_vapour', 'density of saturated water vapour in kg/m³ at 10°C', 'double', 0.0094),
            ('R_co2', 'specific gas constant of CO2 in J/(kg*4.5K)', 'double', 188.91),
            ('lv', 'specific latent evaporation heat(condensation heat) in J/kg', 'double', 2.52e6),
            ('ls', 'specific latent vaporisation heat(sublimation heat) in J/kg', 'double', 2.83e6),
            ('cp_l', 'specific heat capacity of dry air at constant pressure and 20°C in J/(kg K)', 'double', 1005.0),
            ('cv_l', 'specific heat capacity of dry air at constant volume and 20°C in J/(kg K)', 'double', 717.0),
            ('lamda', 'heat transfer coefficient of air in W/(m K)', 'double', 0.0262),
            ('r_co2', 'density of CO2 in kg/m³ at 25°C', 'double', 0.0019767),
            ('gam', 'constant slope of temperature    gam = 6.5 K/1000 m', 'double', 0.0065),

            ('u_0', 'annual mean of surface wind velocity in m/s, 8 m/s compare to 28.8 km/h', 'double', 8.0),
            ('t_00', 'temperature in K compare to -37°C', 'double', 236.15),
            ('t_000', 'temperature in K compare to -20°C', 'double', 253.15),
            ('s_0', 'entropy at 0°C, cp_l * t_0 in J/kg', 'double', 274515.75),
            ('c_0', 'maximum value of water vapour in kg/kg', 'double', 0.035),

#            ('co2_0', 'maximum value of CO2 in ppm at preindustrial times', 'double', 280.0),
            ('co2_0', 'maximum value of CO2 in ppm at preindustrial times', 'double', 380.0),
            ('co2_paleo', 'value at modern times', 'double', 330.0),
            ('co2_tropopause', 'minimum rate CO2 at tropopause 320.0 ppm', 'double', 385.0),
            ('co2_vegetation', 'value compares to ppm of co2 consumed by the vegetation', 'double', 140.0),
            ('co2_ocean', 'value compares to ppm of co2 consumed by the vegetation', 'double', 0.0),
            ('co2_land', 'value compares to ppm of co2 consumed by the vegetation', 'double', 0.0),
            ('co2_scale', 'multiplier applied to the whole CO2 field for sensitivity experiments (1.0 = field as built; 2.0 = doubled CO2)', 'double', 1.0),

            ('c_land', 'water vapour reduction on land(60% of the saturation value)', 'double', 66),
            ('c_ocean', 'water vapour reduction on sea surface(64% of the saturation value)', 'double', 70),

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
