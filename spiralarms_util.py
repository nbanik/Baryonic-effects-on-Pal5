import glob
import pickle
import numpy
import numpy as np
from numpy.polynomial import Polynomial
from scipy import ndimage, signal, interpolate, integrate
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014, turn_physical_off, plotDensities,evaluateDensities,plotPotentials,SpiralArmsPotential
from galpy.util import bovy_conversion, save_pickles, bovy_coords, bovy_plot
from scipy import ndimage, signal, interpolate, integrate,optimize
import pal5_util
import seaborn as sns
import astropy.units as u
from galpy import potential
from matplotlib import cm, pyplot
from galpy.potential import DehnenSmoothWrapperPotential as DehnenWrap
from galpy.potential import SCFPotential
import streamspraydf

ro=8.
vo=220.


def spiral_arms_potential(FR_frac=1.,cos=True,N=2,pat_speed=24.5,pitch_angle=9.9,r_ref=8.,Rs=7.,phi0=26.,H=0.3):
    
    
        phi0=np.radians(phi0)
        omega=pat_speed*(ro/vo)
        alpha=numpy.radians(pitch_angle)
        r_ref/=ro
        Rs/=ro
        H/=ro

        # percentage of the radial force to set the amplitude of the spiral
        FR_frac=FR_frac*0.01  
        
        if cos :
            Cs=[8./(3.*numpy.pi),0.5,8./(15.*numpy.pi)]
            
        else :
            Cs=[1]
            
        #compute max radial force for amp=1
        pp=np.linspace(0.,2.*np.pi,1000)
        FR_nonaxi=[]
        spiral_pot_amp1=SpiralArmsPotential(amp=1.,N=N,omega=omega,alpha=alpha,phi_ref=phi0,r_ref=r_ref,H=H,Rs=Rs,Cs=Cs)
        turn_physical_off(spiral_pot_amp1)
        
        for ii in range(len(pp)):
             FR_nonaxi.append(potential.evaluateRforces(spiral_pot_amp1,R=8./ro,z=0.,phi=pp[ii],t=0.))

        interp_FR_nonaxi= interpolate.InterpolatedUnivariateSpline(pp,FR_nonaxi)

        #fmin, because radial force is negative
        max_phi=optimize.fmin(interp_FR_nonaxi, 0.,disp=0)

        max_FR_nonaxi= interp_FR_nonaxi(max_phi)[0]
                
        #compute the radial force due to the axisymmetric potential
        FR_axi=potential.evaluateRforces(MWPotential2014,R=8./ro,z=0.,t=0.)

        #compute the correct amplitude
        amp= numpy.abs(FR_frac*FR_axi/max_FR_nonaxi)
        
        #setup spiral potential with correct amplitude
        spiralpot=SpiralArmsPotential(amp=amp,N=N,omega=omega,alpha=alpha,phi_ref=phi0,r_ref=r_ref,H=H,Rs=Rs,Cs=Cs)
        turn_physical_off(spiralpot)
        
        MWspiralpot = MWPotential2014 + [spiralpot]
        
        return MWspiralpot


def sample_streamdf_pal5_noprog_spiral(N,spiralpot,nospiralpot=MWPotential2014,fo='spiral_trailing.dat',trailing=True):
    
              
        if trailing :
            sdf_trailing= pal5_util.setup_pal5model(pot=nospiralpot)
            R,vR,vT,z,vz,phi,dt= sdf_trailing.sample(n=N,returndt=True)
            fo=open(fo,'w')
          
        
        else :
            sdf_leading= pal5_util.setup_pal5model(pot=nospiralpot,leading=True)
            R,vR,vT,z,vz,phi,dt= sdf_leading.sample(n=N,returndt=True)
            fo_lead=fo.replace('trailing','leading')
            fo=open(fo_lead,'w')
              
        finalR= numpy.empty(N)
        finalvR=numpy.empty(N)
        finalvT=numpy.empty(N)
        finalvz=numpy.empty(N)
        finalphi= numpy.empty(N)
        finalz= numpy.empty(N)
        tt=numpy.empty(N)

        for ii in range(N):

                o= Orbit([R[ii],vR[ii],vT[ii],z[ii],vz[ii],phi[ii]])
                o.turn_physical_off()
                ts= numpy.linspace(0.,-dt[ii],1001)

                o.integrate(ts,nospiralpot)
                orb=Orbit([o.R(ts[-1]),o.vR(ts[-1]),o.vT(ts[-1]),o.z(ts[-1]),o.vz(ts[-1]),o.phi(ts[-1])])
                                
                ts_future= numpy.linspace(-dt[ii],0.,1001)
                #forward integrate in barred potential
                orb.integrate(ts_future,spiralpot)
                finalR[ii]= orb.R(ts_future[-1])
                finalphi[ii]= orb.phi(ts_future[-1])
                finalz[ii]= orb.z(ts_future[-1])
                finalvR[ii]=orb.vR(ts_future[-1])
                finalvT[ii]=orb.vT(ts_future[-1])
                finalvz[ii]=orb.vz(ts_future[-1])
                tt[ii]=dt[ii]
                
           
        fo.write("#R   phi   z   vR    vT    vz    ts" + "\n")
    
        for jj in range(N):
            fo.write(str(finalR[jj]) + "   " + str(finalphi[jj]) + "   " + str(finalz[jj]) + "   " + str(finalvR[jj]) + "   " + str(finalvT[jj]) + "   " + str(finalvz[jj]) + "   " + str(tt[jj]) + "\n")
        
        fo.close()
    
        return None