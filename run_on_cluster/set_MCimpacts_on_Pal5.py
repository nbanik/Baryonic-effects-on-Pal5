import numpy as np
import pickle
from astropy.io import fits
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from galpy.util import bovy_conversion, bovy_coords, save_pickles, bovy_plot
from galpy.potential import MWPotential2014, turn_physical_off, vcirc
import astropy.units as u
from galpy.orbit import Orbit
from optparse import OptionParser
import GMC_util
import pal5_util


ro=8.
#paper on MC used R0=8.5 kpc, using ro=8. as of now.
vo=220.

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
        
    parser.add_option("--ind",dest='ind',default=None,
                      type='int',
                      help="index of apar")
    return parser
    
parser= get_options()
options,args= parser.parse_args()

timpact,apar,x_stream,y_stream,z_stream,_,_,_=GMC_util.aparxv_stream_from_pkl(sampling=4096,nchunks=64)

impactMC_ind,M_mc,rs_mc,v_mc,impactb,impact_angle,tmin=GMC_util.compute_impact_parameters(timpact,apar,x_stream,y_stream,z_stream,nchunks=64,sampling_low=128,imp_fac=5.,Mmin=10**5.,rand_rotate=True)

#load the lower timpact pkl file
with open('pkl_files/pal5pepper_128sampling_Plummer_MW2014.pkl','rb') as savefile:
            sdf_pepper= pickle.load(savefile,encoding='latin1')
        
sdf_smooth= pal5_util.setup_pal5model()

sdf_pepper.set_impacts(impactb=impactb,subhalovel=v_mc,impact_angle=impact_angle,timpact=tmin,rs=rs_mc,GM=M_mc)

pepperfilename='Pal5_4096_on_128impact_Plummer_Mmin105_MC_rand_rotate_{}.pkl'.format(options.ind)

save_pickles(pepperfilename,sdf_pepper)
