import os, os.path
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
                      
    parser.add_option("--td",dest='td',default=5.,
                      type='float',
                      help="tdisrupt in Gyr")
                      
    parser.add_option("--fo",dest='ft',default=None,
                      help="output file name")
                      
    parser.add_option("--t_on",dest='t_on',default=-5.,
                      type='float',
                      help="time in the past when the bar acquired full strength")
                      
    return parser


ro=8.
vo=220.

def galcencyl_to_lbd(R,phi,Z,degree=True):
    xyz=bovy_coords.galcencyl_to_XYZ(R,phi,Z)
    lbd=bovy_coords.XYZ_to_lbd(xyz[0],xyz[1],xyz[2],degree=degree)
    return lbd[0], lbd[1], lbd[2]
    
parser= get_options()
options,args= parser.parse_args()
    
#pat_speed=[39.,43.,47.]

Mbar=10.**10

ang=27.
td=options.td
t_on=[-4.,-3.,-2.,-1.]
pat_speed=39.

Ac,As=SCFbar_util.compute_Acos_Asin()

N=500


for ii in range(len(t_on)):
    
    
        
    ft=options.ft
    fl=ft.replace('trailing','leading')
    
    folder='/ufrc/tan/nilanjan1/galpy/sampled_SCFbar/1010Msun_39patspeed_{}Gyr_spraydf/trailing/'.format(int(np.abs(t_on[ii])))
    folder_lead=folder.replace('trailing','leading')
    
    ft1=ft.replace('5Gyr',str(int(np.abs(t_on[ii])))+'Gyr')
    fl1=fl.replace('5Gyr',str(int(np.abs(t_on[ii])))+'Gyr')
        
    if not os.path.exists(folder):
        os.makedirs(folder,exist_ok=True)
        
    if not os.path.exists(folder_lead):
        os.makedirs(folder_lead,exist_ok=True)
    
    fpatht=folder + ft1
    fpathl=folder_lead + fl1
    
    barpot,nobarpot=SCFbar_util.Particle_Spray_MWPotentialSCFbar(mbar=Mbar,Acos=Ac,Asin=As,t_on=t_on[ii],pat_speed=pat_speed,tstream=td)
    
    p5= Orbit([229.018,-0.124,23.2,-2.296,-2.257,-58.7],radec=True,ro=ro,vo=vo,solarmotion=[-11.1,24.,7.25])
    
    #convert to galpy units
    pal5=Orbit(p5._orb.vxvv)
    
    #mass of Pal 5 from Dehnen https://arxiv.org/pdf/astro-ph/0401422.pdf
    spdf= streamspraydf.streamspraydf(50000.*units.Msun,progenitor=pal5,pot=barpot,tdisrupt=td*units.Gyr)
    spdft= streamspraydf.streamspraydf(50000.*units.Msun,progenitor=pal5,pot=barpot,leading=False,tdisrupt=td*units.Gyr)
    
    
    
    #Rt,vRt,vTt,zt,vzt,phit
    RvRl,pdtl= spdf.sample(n=N,returndt=True,integrate=True)
    RvRt,pdtt= spdft.sample(n=N,returndt=True,integrate=True)
    
        
    fo_trail=open(fpatht,'w')
    fo_lead=open(fpathl,'w')
    
    fo_trail.write("#R   phi   z   vR    vT    vz    ts" + "\n")
    fo_lead.write("#R   phi   z   vR    vT    vz    ts" + "\n")
    
    for jj in range(N):
            fo_trail.write(str(RvRt[0][jj]) + "   " + str(RvRt[5][jj]) + "   " + str(RvRt[3][jj]) + "   " + str(RvRt[1][jj]) + "   " + str(RvRt[2][jj]) + "   " + str(RvRt[4][jj]) + "   " + str(pdtt[jj]) + "\n")
            fo_lead.write(str(RvRl[0][jj]) + "   " + str(RvRl[5][jj]) + "   " + str(RvRl[3][jj]) + "   " + str(RvRl[1][jj]) + "   " + str(RvRl[2][jj]) + "   " + str(RvRl[4][jj]) + "   " + str(pdtl[jj]) + "\n")
        
    fo_trail.close()
    fo_lead.close()




