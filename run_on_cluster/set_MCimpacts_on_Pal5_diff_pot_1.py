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
import GMC_util_1
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

##########setup Pal5 orbit and potential
paramf=np.genfromtxt(options.param_file,delimiter=',')  

pind=paramf[options.ind][0]

peri=round(paramf[options.ind][1],2)
print (peri)


flat_c=paramf[options.ind][2]
vc=paramf[options.ind][3]
tvo= vc*_REFV0

#indices greater than 14: subtract 1
if pind > 14 :
    pind -=1
    
pind=int(pind)   
potparams_file=np.loadtxt('pot_params.dat',delimiter=',')
potparams=list(potparams_file[pind])
   
pot= MWPotential2014Likelihood.setup_potential(potparams,flat_c,False,False,
                                                       pal5_util_MWfit._REFR0,tvo)


timpact,apar,x_stream,y_stream,z_stream,_,_,_=GMC_util_1.aparxv_stream_from_pkl(pot=pot,sampling=4096,nchunks=64)

impactMC_ind,M_mc,rs_mc,v_mc,impactb,impact_angle,tmin=GMC_util_1.compute_impact_parameters(timpact,apar,x_stream,y_stream,z_stream,pot=pot,nchunks=64,sampling_low=128,imp_fac=5.,Mmin=10**5.,rand_rotate=True)

#load the lower timpact pkl file
with open('pkl_files/pal5pepper_Plummer_128sampling_pot{}_peri{}.pkl'.format(pind,peri),'rb') as savefile:
            sdf_pepper= pickle.load(savefile,encoding='latin1')
        
sdf_smooth= pal5_util.setup_pal5model(pot=pot)

sdf_pepper.set_impacts(impactb=impactb,subhalovel=v_mc,impact_angle=impact_angle,timpact=tmin,rs=rs_mc,GM=M_mc)

pepperfilename='Pal5_4096_on_128impact_Plummer_pot{}_peri{}_Mmin105_MC_rand_rotate_{}.pkl'.format(pind,peri)

save_pickles(pepperfilename,sdf_pepper)
