import matplotlib
matplotlib.use('agg')
import numpy as np
import pickle
import numpy
from astropy.io import fits
from galpy.util import bovy_conversion, bovy_coords, save_pickles, bovy_plot
from galpy.potential import MWPotential2014, turn_physical_off, vcirc
import astropy.units as u
from galpy.orbit import Orbit
import random
import pal5_util_MWfit
import MWPotential2014Likelihood
import os, os.path
import re
import glob
import pickle
import csv
from optparse import OptionParser
_REFR0, _REFV0= MWPotential2014Likelihood._REFR0, MWPotential2014Likelihood._REFV0
ro, vo= _REFR0, _REFV0


def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
    parser.add_option("--ind",dest='ind',default=None,
                      type='int',
                      help="index of potential")
    return parser


def determine_nburn(filename='../pal5_mcmc/mwpot14-fitsigma-0.dat',
                    threshold=0.1,skip=50,
                    return_nsamples=False):
    """Function to detemrine an appropriate nburn for a given chain"""
    # Load the data
    data= numpy.loadtxt(filename,comments='#',delimiter=',')
    lndata= numpy.reshape(data[:,-1],(len(data[:,5])//nwalkers,nwalkers))
    # Perform a running diff wrt skip less
    diff= (lndata-numpy.roll(lndata,skip,axis=0))
    diff[:skip]= -100. # Make sure it's not within the first hundred
    maxln= numpy.nanmax(lndata)
    try:
        indx= (numpy.fabs(numpy.median(diff,axis=1)) < threshold)\
                       *((maxln-numpy.nanmax(lndata,axis=1)) < 1.25)
        if maxln > -22.5:
            indx*= numpy.std(lndata,axis=1) < 3.
        if return_nsamples:
            return len(data)-numpy.arange(len(lndata))[indx][0]*nwalkers
        else:
            return numpy.arange(len(lndata))[indx][0]*nwalkers
    except IndexError:
        if return_nsamples: return 100.
        else: return numpy.prod(lndata.shape)-100


nwalkers= 12

#from each MCMC chain file, pick nsamples
nsamples= 2000

pot_ind=np.arange(0,32,1)
pot_ind=np.delete(pot_ind,14)

t_age= np.linspace(0.,5.,1001)/bovy_conversion.time_in_Gyr(vo,ro)

peri_all=[]

parser= get_options()
options,args= parser.parse_args()

pindx=pot_ind[options.ind]

csvfo= open('pal5_mcmc_selected_chains_pot{}.dat'.format(pindx),'w')
fowriter= csv.writer(csvfo,delimiter=',')

# Load this potential
fn= 'mwpot14-fitsigma-%i.dat' % pindx
with open(fn,'rb') as savefile:
    line1= savefile.readline()
potparams= [float(s) for s in (line1.split(':'.encode())[1].split(','.encode()))]

tnburn= determine_nburn(fn)
tdata= numpy.loadtxt(fn,comments='#',delimiter=',')
tdata= tdata[tnburn::]

rand_indx=random.sample(range(len(tdata)),nsamples)

peri=[]

for jj in rand_indx:
    
    tvo= tdata[jj][1]*_REFV0
    pot= MWPotential2014Likelihood.setup_potential(potparams,tdata[jj][0],False,False,
                                                    pal5_util_MWfit._REFR0,tvo)

    # Now compute the stream model for this setup
    dist= tdata[jj][2]*22.
    pmra= -2.296+tdata[jj][3]+tdata[jj][4]
    pmdecpar= 2.257/2.296
    pmdecperp= -2.296/2.257
    pmdec= -2.257+tdata[jj][3]*pmdecpar+tdata[jj][4]*pmdecperp
    vlos= -58.7
    sigv= 0.4*numpy.exp(tdata[jj][5])
    
    
    prog= Orbit([229.018,-0.124,dist,pmra,pmdec,vlos],
                radec=True,ro=ro,vo=tvo,
                solarmotion=[-11.1,24.,7.25]).flip()
    
    prog.integrate(t_age,pot)
    peri=prog.rperi()
    
    out=[peri,tdata[jj][3],tdata[jj][0],tdata[jj][1],sigv]
    out.extend([229.018,-0.124,dist,pmra,pmdec,vlos])
    
    fowriter.writerow(out)
    
csvfo.flush()
csvfo.close()
