import os, os.path
import glob
import pickle
import numpy
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial
from scipy import ndimage, signal, interpolate, integrate
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,turn_physical_off, MiyamotoNagaiPotential
from galpy.util import bovy_conversion, save_pickles, bovy_coords, bovy_plot
import pal5_util
import seaborn as sns
import astropy.units as u
from galpy import potential
from optparse import OptionParser
import argparse
from galpy.potential import DehnenBarPotential
from galpy.potential import DehnenSmoothWrapperPotential as DehnenWrap
import spiralarms_util
from optparse import OptionParser

def galcencyl_to_lbd(R,phi,Z,degree=True):
    xyz=bovy_coords.galcencyl_to_XYZ(R,phi,Z)
    lbd=bovy_coords.XYZ_to_lbd(xyz[0],xyz[1],xyz[2],degree=degree)
    return lbd[0], lbd[1], lbd[2]

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
        
    parser.add_option("--pat_speed",dest='pat_speed',default=40.,
                      type='float',
                      help="pattern speed of the bar")
                                             
    parser.add_option("--ind",dest='ind',default=None,
                      help="index")
                      
    parser.add_option("--t_on",dest='t_on',default=-5.,
                      type='float',
                      help="time in the past when the bar acquired full strength")
                      
    return parser
    
parser= get_options()
options,args= parser.parse_args()

ro,vo= 8., 220.


Nsamples=500

#pat_speed=[39.,43.,47.]

ind=options.ind
t_on=[-5.,-3.,-1.]

FR_frac=[0.5,1.0]
N=[2,4]
pat_speed=[19.5,24.5]


for t in t_on:
    
    for fr in FR_frac :
        
        for n in N :
            
            for p in pat_speed :
                
        
                    folder='/ufrc/tan/nilanjan1/galpy/sampled_spiral/spiralN{}_{}patspeed_FR{}_{}Gyr_spraydf/trailing/'.format(n,p,fr,int(np.abs(t)))
                    folder_lead=folder.replace('trailing','leading')
                    
                    fot='sample500_spiral_spraydf_trailing_N{}_{}patspeed_{}Gyr_{}FR_{}.dat'.format(n,p,int(np.abs(t)),fr,ind)
                    
                    if not os.path.exists(folder):
                        os.makedirs(folder,exist_ok=True)
                        
                    if not os.path.exists(folder_lead):
                        os.makedirs(folder_lead,exist_ok=True)
                    
                    fpath=folder + fot
                    
                    spiralpot=spiralarms_util.spiral_arms_potential(FR_frac=fr,t_on=t,N=n,pat_speed=p)
                    
                    
                    spiralarms_util.sample_spraydf_pal5_spiral(Nsamples,spiralpot,fo=fpath)
                    spiralarms_util.sample_spraydf_pal5_spiral(Nsamples,spiralpot,trailing=False,fo=fpath)
                    
                    







