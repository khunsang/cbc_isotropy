from optparse import OptionParser
#from ligo.lw import ligolw
#from ligo.lw import lsctables
from glue.lal import Cache
#from ligo.lw import utils as ligolw_utils
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ligolw, table, lsctables
from astropy.units.si import sday
from bilby.gw import conversion
import os
import sys
import numpy as np
import lalsimulation as lalsim
import pandas as pd
import lal

@lsctables.use_in
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass

def gmst_rad_from_geocetric_end_time(gps_time):
    phase = (float(gps_time) / sday.si.scale ) * 2 * np.pi
    gmst_rad = phase % (2*np.pi)
    return gmst_rad

def logitude_to_RA_no_lalsuite(longitude, gps_time):
    gmst = gmst_rad_from_geocetric_end_time(gps_time)
    ra = longitude + gmst
    return ra % (2*np.pi)


def parse_command_line():
    parser = OptionParser(description = __doc__)
    parser.add_option("--injection-file", metavar = "filename", default= "full_injections.xml.gz", help = "Set the injection xml file.")
    parser.add_option("--eos", metavar = "eos name", default="APR", help = "EOS for BNS params.")
    parser.add_option("--flow", metavar = "flow Hz", default="9", help = "lower cut-off frequency", type='float')
    parser.add_option("--inj-dir", metavar = "injection directory", default=".", help = "Set the injection directory")
    options, filenames = parser.parse_args()
    return options, filenames





options, filenames = parse_command_line()

inj_dir = options.inj_dir


inj_xml_file = os.path.join( inj_dir, options.injection_file)



cmd_injection_file = 'lalapps_inspinj --m-distr componentMass --min-mass1 5.0 --max-mass1 65.0 \
        --min-mass2 3.0 --max-mass2 55.0 --min-mtotal 8.0 --max-mtotal 120.0\
        --enable-spin --aligned  --min-spin1 0  --min-spin2 0  --max-spin1 0.9 --max-spin2 0.9\
        --i-distr uniform --l-distr random --d-distr volume --min-distance 40000 --max-distance 2500000\
        --seed 100 --t-distr uniform --time-step 2000 --time-interval 25 --gps-start-time 1126051217  --gps-end-time 1128299417\
        --f-lower {} --waveform TaylorF2 --output {}'.format( options.flow, inj_xml_file )


os.system(cmd_injection_file)


#Read inj files sim_inspiral
xmldoc = ligolw_utils.load_filename( inj_xml_file, verbose = True, contenthandler = LIGOLWContentHandler)

xml_table = table.get_table(  xmldoc, lsctables.SimInspiralTable.tableName  )


mass1_samples = xml_table.get_column('mass1')
mass2_samples = xml_table.get_column('mass2')

mchirp_samples = xml_table.get_column('mchirp')
mass_ratio_samples = conversion.component_masses_to_mass_ratio(mass1_samples, mass2_samples)

geocent_end_time_samples =  xml_table.get_column('geocent_end_time')


longitude_samples = xml_table.get_column('longitude')
### make sure that longitudes are in the range [0,2\pi]

ra_samples = []

for i in range(len(geocent_end_time_samples)):
    ra = logitude_to_RA_no_lalsuite(  longitude_samples[i],
            geocent_end_time_samples[i] )
    ra_samples.append(ra)

ra_samples = np.asarray(ra_samples)
dec_samples = xml_table.get_column('latitude')
distance_samples = xml_table.get_column('distance')


spin1z_samples = xml_table.get_column('spin1z')
spin2z_samples = xml_table.get_column('spin2z')



inj_dict = { "mass_1":mass1_samples, "mass_2":mass2_samples,  
        "luminosity_distance": distance_samples, "dec": dec_samples, "ra": ra_samples,
        "psi":xml_table.get_column('polarization'), "phase": xml_table.get_column('coa_phase'),
        "theta_jn": xml_table.get_column('inclination'),  "geocent_time":geocent_end_time_samples}


injection_set = pd.DataFrame.from_dict(inj_dict)
injection_set.to_csv( os.path.join( inj_dir, "injection_cbc_params.dat"), sep=' ', index=False)
