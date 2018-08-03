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
import GMC_util
import MWPotential2014Likelihood
_REFR0, _REFV0= MWPotential2014Likelihood._REFR0, MWPotential2014Likelihood._REFV0
ro, vo= _REFR0, _REFV0



def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
                          
    parser.add_option("--chain_ind",dest='chain_ind',default=None,
                      type='int',
                      help="index of line from parama_file to be used to compute streampeppered object")
    
    #set timpact chunks                        
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

            

#sdf_smooth= pal5_util.setup_pal5model()
pepperfilename= 'pkl_files/pal5pepper_Plummer_{}sampling_chainind{}_{}.pkl'.format(len(timpacts),options.chain_ind,options.tind)

sdf_pepper=GMC_util.make_nondefault_pal5stream(options.chain_ind,timpact=timpactn)

save_pickles(pepperfilename,sdf_pepper)
