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


def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
    parser.add_option("--ind",dest='ind',default=None,
                      type='int',
                      help="index of apar")
    return parser
    


            
apar=numpy.arange(0.,1.75,0.01)

parser= get_options()
options,args= parser.parse_args()

for ii in range(1,21):
    with open('Pal5_4096_on_128impact_Plummer_Mmin105_MC_rand_rotate_{}.pkl'.format(ii),'rb') as savefile:
            #sdf_smooth=pickle.load(savefile,encoding='latin1')
            sdf_pepper= pickle.load(savefile,encoding='latin1')

    sdf_smooth= pal5_util.setup_pal5model()
    a=apar[options.ind]
    dens_unp= sdf_smooth._density_par(a)
    omega_unp=sdf_smooth.meanOmega(a,oned=True)
    
    densOmega= sdf_pepper._densityAndOmega_par_approx(a)


    fo=open('dens_Omega/densOmega_4096_on_128_Plummer_Mmin105_rand_rotate{}_{}.dat'.format(ii,options.ind),'w')
    fo.write('#apar   dens_unp   dens  omega_unp   omega' + '\n')
    fo.write(str(a) + '   ' + str(dens_unp) + '   ' + str(densOmega[0]) + '   ' + str(omega_unp) + '   ' + str(densOmega[1]) )
    
    fo.close()

