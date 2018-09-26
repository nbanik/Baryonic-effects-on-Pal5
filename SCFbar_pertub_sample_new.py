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
from gd1_util import R0, V0
import custom_stripping_df
import seaborn as sns
import astropy.units as u
from galpy import potential
from optparse import OptionParser
import argparse
from galpy.potential import DehnenBarPotential
from galpy.potential import DehnenSmoothWrapperPotential as DehnenWrap
import SCFbar_util

def galcencyl_to_lbd(R,phi,Z,degree=True):
    xyz=bovy_coords.galcencyl_to_XYZ(R,phi,Z)
    lbd=bovy_coords.XYZ_to_lbd(xyz[0],xyz[1],xyz[2],degree=degree)
    return lbd[0], lbd[1], lbd[2]

parser = argparse.ArgumentParser(description='My app description')
parser.add_argument('-o', '--output', help='Path to output file')
args = parser.parse_args()

fo=args.output

ro,vo= 8., 220.

Ac,As=SCFbar_util.compute_Acos_Asin()

Mbar=10**10.
pat_speed=40.
ang=27.

barpot,nobarpot=SCFbar_util.Particle_Spray_MWPotentialSCFbar(mbar=Mbar,Acos=Ac,Asin=As,t_on=-2.)

SCFbar_util.sample_perturbed_pal5_new(100,barpot,nobarpot,fo=fo)
SCFbar_util.sample_perturbed_pal5_new(100,barpot,nobarpot,trailing=False,fo=fo)





