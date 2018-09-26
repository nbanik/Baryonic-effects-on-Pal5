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

Mbar=[6.*10**9.,8*10**9.,1.2*10**10.,1.4*10**10]

Mbar_str=['6_10_9','8_10_9','12_10_9','14_10_9']

pat_speed=options.pat_speed
ft=options.ft
t_on=options.t_on
ang=27.

for ii in range(len(Mbar)):
    
    folder='/ufrc/tan/nilanjan1/galpy/sampled_SCFbar/' + Mbar_str[ii] + 'Msun_{}_nobarpot_streamdf/trailing/'.format(int(pat_speed))
    folder_lead=folder.replace('trailing','leading')
    ft1=ft.replace('1010',Mbar_str[ii])
    if not os.path.exists(folder):
        os.makedirs(folder)
        
    if not os.path.exists(folder_lead):
        os.makedirs(folder_lead)
    
    fpath=folder + ft1
    
    barpot,nobarpot=SCFbar_util.Particle_Spray_MWPotentialSCFbar(mbar=Mbar[ii],Acos=Ac,Asin=As,t_on=t_on,pat_speed=pat_speed)


    SCFbar_util.sample_streamdf_pal5_noprog(N,nobarpot,nobarpot,fo=fpath)
    SCFbar_util.sample_streamdf_pal5_noprog(N,nobarpot,nobarpot,trailing=False,fo=fpath)








