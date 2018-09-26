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
import SCFbar_util
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
                                             
    parser.add_option("--fo",dest='ft',default=None,
                      help="trailing file name")
                      
    parser.add_option("--t_on",dest='t_on',default=-5.,
                      type='float',
                      help="time in the past when the bar acquired full strength")
                      
    return parser
    
parser= get_options()
options,args= parser.parse_args()

ro,vo= 8., 220.

Ac,As=SCFbar_util.compute_Acos_Asin()

N=500

#pat_speed=[39.,43.,47.]

Mbar=10.**10

ft=options.ft
t_on=[-4.,-3.,-2.,-1.]

pat_speed=39.
ang=27.


for ii in range(len(t_on)):
        
    folder='/ufrc/tan/nilanjan1/galpy/sampled_SCFbar/1010Msun_39patspeed_{}Gyr_streamdf/trailing/'.format(int(np.abs(t_on[ii])))
    folder_lead=folder.replace('trailing','leading')
    
    ft1=ft.replace('5Gyr',str(int(np.abs(t_on[ii])))+'Gyr')
    
    if not os.path.exists(folder):
        os.makedirs(folder,exist_ok=True)
        
    if not os.path.exists(folder_lead):
        os.makedirs(folder_lead,exist_ok=True)
    
    fpath=folder + ft1
    
    barpot,nobarpot=SCFbar_util.Particle_Spray_MWPotentialSCFbar(mbar=Mbar,Acos=Ac,Asin=As,t_on=t_on[ii],pat_speed=pat_speed)
    
    
    SCFbar_util.sample_streamdf_pal5_noprog(N,barpot,nobarpot,fo=fpath)
    SCFbar_util.sample_streamdf_pal5_noprog(N,barpot,nobarpot,trailing=False,fo=fpath)








