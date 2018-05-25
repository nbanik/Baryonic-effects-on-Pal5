import numpy
from scipy import ndimage, signal, interpolate

from scipy import integrate
#interpolate to get dens(apar)

dat=numpy.loadtxt('apars_dens_simind25.dat')
apars=dat[:,0]
dens=dat[:,1]

ipdens= interpolate.InterpolatedUnivariateSpline(apars,dens)

norm=integrate.romberg(ipdens,min(apars),max(apars),divmax=20)

def dens_pdf(a):
    
    return ipdens(a)/norm

def dens_cdf(a):
      
    return integrate.romberg(dens_pdf,min(apars),a,divmax=20)
    
from scipy.optimize import brentq
def icdf(a):
        return dens_cdf(a) - numpy.random.uniform()

N=1000
sample_apar= numpy.empty(N)

for ii in range(N):
    sample_apar[ii]=brentq(icdf,min(apars),max(apars))
    
    
fo=open('sample_apars_2.dat','w')

for i in range(N):
    fo.write(str(sample_apar[i]) + '\n')
    
    
fo.close()
