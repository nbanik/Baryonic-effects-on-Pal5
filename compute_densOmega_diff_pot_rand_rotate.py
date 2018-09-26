# script to run simulations of stream peppering
import os, os.path
import csv
import time
import pickle
import numpy
import matplotlib
matplotlib.use('agg')
from scipy import integrate, interpolate
from optparse import OptionParser
from galpy.util import bovy_conversion
import gd1_util
import pal5_util
from gd1_util import R0,V0
from scipy.integrate import quad
from scipy.optimize import brentq
import GMC_util


def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
    parser.add_option("--ind",dest='ind',default=None,
                      type='int',
                      help="index of apar")
                      
    parser.add_option("--chain_ind",dest='chain_ind',default=None,
                      type='int',
                      help="chain_index")
                      
    parser.add_option("--td",dest='td',default=5.,
                      type='float',
                      help="tdisrupt in Gyr")
                      
    parser.add_option("--nsampling",dest='nsam',default=4096,
                      type='int',
                      help="no of samplings in the hi res")
    return parser
    


            
apar=numpy.arange(0.,2.0,0.01)

parser= get_options()
options,args= parser.parse_args()

for ii in range(20):
    
    with open('Pal5_{}_on_128impact_Plummer_td{}_Mmin105_chainind_{}_seedind{}.pkl'.format(options.nsam,options.td,options.chain_ind,ii),'rb') as savefile:
        #sdf_smooth=pickle.load(savefile,encoding='latin1')
        sdf_pepper= pickle.load(savefile,encoding='latin1')

    sdf_smooth= GMC_util.make_nondefault_pal5stream(options.chain_ind,td=options.td)
    
    a=apar[options.ind]
    dens_unp= sdf_smooth._density_par(a)
    omega_unp=sdf_smooth.meanOmega(a,oned=True)
    
    densOmega= sdf_pepper._densityAndOmega_par_approx(a)
    
    
    fo=open('dens_Omega/densOmega_{}_on_128_Plummer_td{}_chainind{}_seedind{}_Mmin105_{}.dat'.format(options.nsam,options.td,options.chain_ind,ii,options.ind),'w')
    fo.write('#apar   dens_unp   dens  omega_unp   omega' + '\n')
    fo.write(str(a) + '   ' + str(dens_unp) + '   ' + str(densOmega[0]) + '   ' + str(omega_unp) + '   ' + str(densOmega[1]) )
    
    fo.close()

