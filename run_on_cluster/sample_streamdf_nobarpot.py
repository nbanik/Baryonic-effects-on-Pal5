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
                                             
    parser.add_option("--fo",dest='fo',default=None,
                      help="trailing file name")
                      
    parser.add_option("--t_on",dest='t_on',default=-5.,
                      type='float',
                      help="time in the past when the bar acquired full strength")
                      
    return parser
    
parser= get_options()
options,args= parser.parse_args()

ro,vo= 8., 220.

Ac,As=SCFbar_util.compute_Acos_Asin()

Mbar=10**10.
pat_speed=40.
fot=options.fo
t_on=options.t_on
ang=27.

barpot,nobarpot=SCFbar_util.Particle_Spray_MWPotentialSCFbar(mbar=Mbar,Acos=Ac,Asin=As,t_on=t_on,pat_speed=pat_speed)



def sample_streamdf_pal5_nobarpot(nobarpot,fot,N=500,trailing=True):

        if trailing :

            sdf_trailing= pal5_util.setup_pal5model(pot=nobarpot)

            R,vR,vT,z,vz,phi,dt= sdf_trailing.sample(n=N,returndt=True)
            fo='/ufrc/tan/nilanjan1/galpy/sampled_SCFbar/testpeak_1010Msun_nobarpot_streamdf/trailing/' + fot
            fo=open(fo,'w')

        else :

            sdf_leading= pal5_util.setup_pal5model(pot=nobarpot,leading=True)

            R,vR,vT,z,vz,phi,dt= sdf_leading.sample(n=N,returndt=True)

            fol='/ufrc/tan/nilanjan1/galpy/sampled_SCFbar/testpeak_1010Msun_nobarpot_streamdf/trailing/' + fot
            fo_lead=fol.replace('trailing','leading')

            fo=open(fo_lead,'w')

                        

        fo.write("#R   phi   z   vR    vT    vz    ts" + "\n")
        for jj in range(N):

            fo.write(str(R[jj]) + "   " + str(phi[jj]) + "   " + str(z[jj]) + "   " + str(vR[jj]) + "   " + str(vT[jj]) + "   " + str(vz[jj]) + "   " + str(dt[jj]) + "\n")

        

        fo.close()

    

        return None



sample_streamdf_pal5_nobarpot(nobarpot,fot=fot,N=500,trailing=True)
sample_streamdf_pal5_nobarpot(nobarpot,fot=fot,N=500,trailing=False)




