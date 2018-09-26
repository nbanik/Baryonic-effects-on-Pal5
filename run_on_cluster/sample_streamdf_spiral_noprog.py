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
                      help="pattern speed of the spiral")
                                             
    parser.add_option("--fo",dest='fo',default=None,
                      help="trailing file name")
                      
    parser.add_option("-N",dest='N',default=2,
                      type='int',
                      help="no of arms")
                      
    parser.add_option("--pitch_angle",dest='pitch_angle',default=9.9,
                      type='float',
                      help="pitch angle of the spiral arms")
                      
    parser.add_option("--FR_frac",dest='FR_frac',default=1.,
                      type='float',
                      help="percentage of the radial force")
                      
    return parser
    
parser= get_options()
options,args= parser.parse_args()

ro,vo= 8., 220.

FR_frac=options.FR_frac
N=options.N
pat_speed=options.pat_speed
fo=options.fo
pitch_angle=options.pitch_angle

nsample=1000

MWspiralpot=spiralarms_util.spiral_arms_potential(FR_frac=FR_frac,cos=True,N=N,pat_speed=pat_speed,pitch_angle=pitch_angle,r_ref=8.,Rs=7.,phi0=26.,H=0.3)

spiralarms_util.sample_streamdf_pal5_noprog_spiral(nsample,spiralpot=MWspiralpot,nospiralpot=MWPotential2014,fo=fo)
spiralarms_util.sample_streamdf_pal5_noprog_spiral(nsample,spiralpot=MWspiralpot,nospiralpot=MWPotential2014,fo=fo,trailing=False)












