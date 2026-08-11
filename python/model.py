#!/usr/bin/env python

from pyatom import Atmosphere, Hydrosphere


class Model(object):
    """
    ATOM Model Object
    """
    def __init__(self):
        self.atm = Atmosphere()
        self.hyd = Hydrosphere()
        self.time_slice = ""
        self.config_atm_xml = "config_atm_xml"
        self.config_hyd_xml = "config_hyd_xml"


    def print_config_atm(self, config_atm_xml):
        print("\n\n\n\n   Paleo Atmospheric Circulations") 
        print("\n\n   atmosphere configuration file name is         ", self.config_atm_xml)
        print("   output path for config_atm_xml file is        ", self.atm.config_xml_path.decode('utf-8'))
        print("   output path for atmosphere results is         ", self.atm.output_path.decode('utf-8'))
        print("   topography path is at                         ", self.atm.bathymetry_path.decode('utf-8'))

    def print_config_hyd(self, config_hyd_xml):
        print("\n\n\n\n   Paleo Ocean Circulations") 
        print("\n\n   ocean configuration file name is              ", self.config_hyd_xml)
        print("   output path for config_atm_xml file is at     ", self.hyd.config_xml_path.decode('utf-8'))
        print("   output path for ocean results is at           ", self.hyd.output_path.decode('utf-8'))
        print("   bathymetry path is at                         ", self.hyd.bathymetry_path.decode('utf-8'))

    def run_Model_atm(self, t_s):
        self.time_slice = t_s
        print("\n   run_Model for the Paleo Atmospheric Circulations code prepared for time-slice    Ma = ", self.time_slice)
        self.atm.run()
        print("\n    successfully terminated Paleo Atmospheric Circulations code for time-slice    Ma = ", self.time_slice)
        print("\n")

    def run_Model_hyd(self, t_s):
        self.time_slice = t_s
        print("\n   run_Model for the Paleo Ocean Circulations code prepared for time-slice    Ma = ", self.time_slice)
        self.hyd.run()
        print("\n    successfully terminated Paleo Ocean Circulations code for time-slice    Ma = ", self.time_slice)
        print("\n")


#atm = Model()
#atm.print_config_atm("config_atm_xml")
#atm.run_Model_atm(0)


hyd = Model()
hyd.print_config_hyd("config_hyd_xml")
hyd.run_Model_hyd(0)
