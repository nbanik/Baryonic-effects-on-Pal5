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
import pal5_util_MWfit
import MWPotential2014Likelihood
_REFR0, _REFV0= MWPotential2014Likelihood._REFR0, MWPotential2014Likelihood._REFV0

ro=8.
#paper on MC used R0=8.5 kpc, using ro=8. as of now.
vo=220.

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
        
    parser.add_option("--ind",dest='ind',default=None,
                      type='int',
                      help="index of line from parama_file to be used to compute streampeppered object")
                      
    parser.add_option("--param_file",dest='param_file',default='rperi_grid_select.dat',
                      help="name of selected parameters")
    return parser
    
parser= get_options()
options,args= parser.parse_args()

sample_low='pkl_files/pal5pepper_Plummer_128sampling_chainind{}.pkl'.format(options.ind)

timpact,apar,x_stream,y_stream,z_stream,vx_stream,vy_stream,vz_stream=GMC_util.aparxv_stream_from_multiple_pkl(pot=options.ind,sampling=4096,npart=64)
impactMC_ind,M_mc,rs_mc,v_mc,impactb,impact_angle,tmin=GMC_util.compute_impact_parameters_GMC(timpact,apar,x_stream,y_stream,z_stream,pot=options.ind,sampling_low_file=sample_low,Mmin=10**5.)

#load the lower timpact pkl file
with open(sample_low,'rb') as savefile:
        sdf_pepper= pickle.load(savefile,encoding='latin1')
        
sdf_pepper.set_impacts(impactb=impactb,subhalovel=v_mc,impact_angle=impact_angle,timpact=tmin,rs=rs_mc,GM=M_mc)

pepperfilename='Pal5_4096_on_128impact_Plummer_Mmin105_chainind_{}.pkl'.format(options.ind)

save_pickles(pepperfilename,sdf_pepper)
