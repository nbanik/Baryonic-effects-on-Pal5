import os, os.path
import pickle
import numpy
import matplotlib
matplotlib.use('agg')
from scipy import interpolate
from galpy.util import bovy_conversion, bovy_plot, save_pickles
import gd1_util
from gd1_util import R0, V0
import seaborn as sns
from matplotlib import cm, pyplot
import simulate_streampepper
from scipy import integrate, interpolate
from optparse import OptionParser
from galpy.util import bovy_conversion
import pal5_util
import pal5_util_MWfit
import MWPotential2014Likelihood
_REFR0, _REFV0= MWPotential2014Likelihood._REFR0, MWPotential2014Likelihood._REFV0
ro, vo= _REFR0, _REFV0



def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
                          
    parser.add_option("--ind",dest='ind',default=None,
                      type='int',
                      help="index of line from parama_file to be used to compute streampeppered object")
                      
    parser.add_option("--param_file",dest='param_file',default='rperi_grid_select.dat',
                      help="name of selected parameters")
                      
    parser.add_option("-t","--timpacts",dest='timpacts',default=None,
                      help="Impact times in Gyr to consider; should be a comma separated list")
    
    parser.add_option("--npart",dest='n',default=4,
                      type='int',
                      help="no of subsets of full timpacts")
                      
    parser.add_option("--tind",dest='tind',default=None,
                      type='int',
                      help="index of timpact chunk")
                      
    
    return parser

    
parser= get_options()
options,args= parser.parse_args()

########setup timpact chunks
def parse_times(times,age):
    if 'sampling' in times:
        nsam= int(times.split('sampling')[0])
        return [float(ti)/bovy_conversion.time_in_Gyr(V0,R0)
                for ti in numpy.arange(1,nsam+1)/(nsam+1.)*age]
    return [float(ti)/bovy_conversion.time_in_Gyr(V0,R0)
            for ti in times.split(',')]
            
timpacts= parse_times(options.timpacts,5.)

nchunk= options.n

ntim=int(len(timpacts)/nchunk)

tim_array=[]

for ii in range(nchunk):
    tim_array.append(timpacts[ii*ntim:(ii+1)*ntim])
    
timpactn=tim_array[options.tind]

##########setup Pal5 orbit and potential
paramf=numpy.genfromtxt(options.param_file,delimiter=',')  

pind=paramf[options.ind][0]

peri=round(paramf[options.ind][1],2)
print (peri)


flat_c=paramf[options.ind][2]
vc=paramf[options.ind][3]
orb=list(paramf[options.ind][4:]) 
tvo= vc*_REFV0

#indices greater than 14: subtract 1
if pind > 14 :
    pind -=1
    
pind=int(pind)   
potparams_file=numpy.loadtxt('pot_params.dat',delimiter=',')
potparams=list(potparams_file[pind])
   
pot= MWPotential2014Likelihood.setup_potential(potparams,flat_c,False,False,
                                                       pal5_util_MWfit._REFR0,tvo)

            

#sdf_smooth= pal5_util.setup_pal5model()
pepperfilename= 'pkl_files/pal5pepper_Plummer_{}sampling_pot{}_peri{}_{}.pkl'.format(len(timpacts),pind,peri,options.tind)
sdf_pepper= pal5_util.setup_pal5model(timpact=timpactn,pot=pot,hernquist=False)
save_pickles(pepperfilename,sdf_pepper)
