import matplotlib
matplotlib.use('agg')
import numpy as np
from galpy.potential import LogarithmicHaloPotential
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014, turn_physical_off, MiyamotoNagaiPotential, plotDensities,evaluateDensities
from galpy.util import bovy_conversion, save_pickles, bovy_coords, bovy_plot
from galpy.df import streamdf,streamgapdf  
from galpy.util import bovy_coords, bovy_conversion
from galpy import potential
from matplotlib import cm, pyplot
import SCFbar_util
from astropy import units
import astropy.units as u
import streamspraydf
import argparse
from galpy.potential import DehnenBarPotential
from galpy.potential import DehnenSmoothWrapperPotential as DehnenWrap

def galcencyl_to_lbd(R,phi,Z,degree=True):
    xyz=bovy_coords.galcencyl_to_XYZ(R,phi,Z)
    lbd=bovy_coords.XYZ_to_lbd(xyz[0],xyz[1],xyz[2],degree=degree)
    return lbd[0], lbd[1], lbd[2]

parser = argparse.ArgumentParser(description='My app description')
parser.add_argument('-o', '--output', help='Path to output file')
args = parser.parse_args()


ro=8.
vo=220.

def galcencyl_to_lbd(R,phi,Z,degree=True):
    xyz=bovy_coords.galcencyl_to_XYZ(R,phi,Z)
    lbd=bovy_coords.XYZ_to_lbd(xyz[0],xyz[1],xyz[2],degree=degree)
    return lbd[0], lbd[1], lbd[2]
    
Mbar=10**10.
pat_speed=40.
ang=27.

Ac,As=SCFbar_util.compute_Acos_Asin()
barpot,nobarpot=SCFbar_util.Particle_Spray_MWPotentialSCFbar(mbar=Mbar,Acos=Ac,Asin=As,t_on=-2.)

p5= Orbit([229.018,-0.124,23.2,-2.296,-2.257,-58.7],radec=True,ro=ro,vo=vo,solarmotion=[-11.1,24.,7.25])

#convert to galpy units
pal5=Orbit(p5._orb.vxvv)

#mass of Pal 5 from Dehnen https://arxiv.org/pdf/astro-ph/0401422.pdf
spdf= streamspraydf.streamspraydf(60000.*units.Msun,progenitor=pal5,pot=nobarpot,tdisrupt=5.*units.Gyr)
spdft= streamspraydf.streamspraydf(60000.*units.Msun,progenitor=pal5,pot=nobarpot,leading=False,tdisrupt=5.*units.Gyr)

N=1000

#Rt,vRt,vTt,zt,vzt,phit
RvRl,pdtl= spdf.sample(n=N,returndt=True,integrate=True)
RvRt,pdtt= spdft.sample(n=N,returndt=True,integrate=True)


ftrail=args.output
flead=ftrail.replace('trailing','leading')

fo_trail=open(ftrail,'w')
fo_lead=open(flead,'w')

fo_trail.write("#R   phi   z   vR    vT    vz    ts" + "\n")
fo_lead.write("#R   phi   z   vR    vT    vz    ts" + "\n")

for jj in range(N):
        fo_trail.write(str(RvRt[0][jj]) + "   " + str(RvRt[5][jj]) + "   " + str(RvRt[3][jj]) + "   " + str(RvRt[1][jj]) + "   " + str(RvRt[2][jj]) + "   " + str(RvRt[4][jj]) + "   " + str(pdtt[jj]) + "\n")
        fo_lead.write(str(RvRl[0][jj]) + "   " + str(RvRl[5][jj]) + "   " + str(RvRl[3][jj]) + "   " + str(RvRl[1][jj]) + "   " + str(RvRl[2][jj]) + "   " + str(RvRl[4][jj]) + "   " + str(pdtl[jj]) + "\n")
    
fo_trail.close()
fo_lead.close()





