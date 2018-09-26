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
    
        
    parser.add_option("--pat_speed",dest='pat_speed',default=24.5,
                      type='float',
                      help="pattern speed of the bar")
                                             
    parser.add_option("--fo",dest='fo',default=None,
                      help="trailing file name")
                      
    parser.add_option("--t_on",dest='t_on',default=-5.,
                      type='float',
                      help="time in the past when the bar acquired full strength")
                      
    parser.add_option("--FR_frac",dest='FR_frac',default=1.,
                      type='float',
                      help="percentage of radial force")
                      
    parser.add_option("-N",dest='N',default=4,
                      type='int',
                      help="No of arms")
                      
    return parser
    
parser= get_options()
options,args= parser.parse_args()

pat_speed=options.pat_speed
fo=options.fo
t_on=options.t_on
FR_frac=options.FR_frac
N=options.N

spiralpot=spiralarms_util.spiral_arms_potential(N=N,FR_frac=FR_frac,t_on=t_on,pat_speed=pat_speed)

Nsamples=500


folder='/ufrc/tan/nilanjan1/galpy/sampled_spiral/spiralN{}_{}patspeed_FR{}_streamdf/trailing/'.format(N,int(pat_speed),FR_frac)
folder_lead=folder.replace('trailing','leading')


if not os.path.exists(folder):
    os.makedirs(folder,exist_ok=True)
    
if not os.path.exists(folder_lead):
    os.makedirs(folder_lead,exist_ok=True)
    
fpath=folder + fo


spiralarms_util.sample_streamdf_pal5_noprog_spiral(Nsamples,spiralpot=spiralpot,fo=fpath)
spiralarms_util.sample_streamdf_pal5_noprog_spiral(Nsamples,spiralpot=spiralpot,trailing=False,fo=fpath)








