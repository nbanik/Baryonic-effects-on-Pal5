import numpy as np
import pickle
from astropy.io import fits
import matplotlib.pyplot as plt
from galpy.util import bovy_conversion, bovy_coords, save_pickles, bovy_plot
from galpy.potential import MWPotential2014, turn_physical_off, vcirc
import astropy.units as u
from galpy.orbit import Orbit
import pal5_util

ro=8.
#paper on MC used R0=8.5 kpc, using ro=8. as of now.
vo=220.

def lbd_to_galcencyl(l,b,d,degree=True):
    xyz=bovy_coords.lbd_to_XYZ(l,b,d,degree=degree)
    Rphiz=bovy_coords.XYZ_to_galcencyl(xyz[:,0],xyz[:,1],xyz[:,2],Xsun=1.,Zsun=0.)
    
    return (Rphiz[:,0],Rphiz[:,1],Rphiz[:,2])

def add_MCs(Mmin=10**6.,rand_rotate=False,vo=vo,ro=ro):
    
    '''
    Setup Molecular clouds and return their
    M,rs,R,vR,vT,z,vz,phi today
    rand_rotate : if True, add random rotation to phi between [0,2pi]
    '''
    
    hdulist=fits.open('molecular_clouds/J_ApJ_834_57_table1.dat.gz.fits')
    aa=hdulist[1].data
    l=aa['GLON']
    b=aa['GLAT']
    #Near or far distance flag (0=near; 1=far) 
    flag=aa['INF']
    Dnear=aa['Dnear']
    Dfar=aa['Dfar']
    znear=aa['znear']
    zfar=aa['zfar']
    #R_sph_gal=aa['Rgal']
    Rnear=aa['Rnear']
    Rfar=aa['Rfar']
    Mnear=aa['Mnear']
    Mfar=aa['Mfar']
    
    D_all=np.empty(len(l))
    zfile=np.empty(len(l))
    rs_all=np.empty(len(l))
    M_all=np.empty(len(l))


    for ii in range(len(l)):
        if flag[ii] == 0 :
            D_all[ii]=Dnear[ii]
            zfile[ii]=znear[ii]
            #rs_all[ii]=Rnear[ii]*0.001  #convert to kpc
            M_all[ii]=Mnear[ii]
        
        
        else :
            D_all[ii]=Dfar[ii]
            zfile[ii]=zfar[ii]
            #rs_all[ii]=Rfar[ii]*0.001 #convert to kpc
            M_all[ii]=Mfar[ii]
    

    R_all,phi_all,z_all= lbd_to_galcencyl(l,b,D_all*(8.5/8.))
    
    def rs(M):
        return 0.1*(M/10**7.)**0.5

    for ii in range(len(M_all)):
        rs_all[ii]=rs(M_all[ii])


    R_all/=ro
    z_all/=ro
    
    M=[]
    rs=[]
    z=[]
    R=[]
    phi=[]

    for ii in range(len(l)):
        if M_all[ii] >= Mmin :
            M.append(M_all[ii])
            rs.append(rs_all[ii])
            z.append(z_all[ii])
            phi.append(phi_all[ii])
            R.append(R_all[ii])
        
    M=np.array(M)
    rs=np.array(rs)
    z=np.array(z)
    R=np.array(R)
    phi=np.array(phi)
    
    if rand_rotate :
        phi+=2*np.pi*np.random.uniform(low=0.,high=1.,size=len(phi))
    
    vT=np.empty(len(M))

    for ii in range(len(M)):
        vT[ii]=vcirc(MWPotential2014,R[ii])
        
    vR=np.zeros(len(M))
    vz=np.zeros(len(M))
    
    coord=[]
    
    for ii in range(len(M)):
        coord.append([R[ii],vR[ii],vT[ii],z[ii],vz[ii],phi[ii]])
    
    return (M,rs,coord)
    
    
    
def compute_min_separation(x_mc,y_mc,z_mc,apar,x_stream,y_stream,z_stream):
    '''
    given (x,y,z) of each molecular cloud, compute the minimum separation from the stream chunks
    
    input: x_mc,y_mc,z_mc of the MCs,
    x_stream,y_stream,z_stream as arrays of the stream chunks
    apar of the stream chunks, in order to output the apar at which the minimum separation from the 
    MC occured
    '''
    
    diffx=x_stream - x_mc
    diffy=y_stream - y_mc
    diffz=z_stream - z_mc
    
    diffxyz=np.c_[diffx,diffy,diffz]
    
    norm = np.linalg.norm(diffxyz,axis=1)
    
    #print (diffx)
    
    #print (norm)
    
    #print (len(x_stream), len(norm))
    
    min_ind=np.argmin(norm)
    
    min_sep= norm[min_ind]
    
    apar_min=apar[min_ind]
    
    return (min_sep,apar_min)
    

def aparxv_stream(sdf_smooth,sdf_pepper):
    
    timpact=sdf_pepper._timpact
    
    apar_full=[]
    x_full=[]
    y_full=[]
    z_full=[]
    vx_full=[]
    vy_full=[]
    vz_full=[]
    
    for kk in range(len(timpact)):
    
        apar=[]
        x=[]
        y=[]
        z=[]
        vx=[]
        vy=[]
        vz=[]
    
        a= sdf_pepper._sgapdfs_coordtransform[timpact[kk]]._kick_interpolatedObsTrackXY
        apar_all=sdf_pepper._sgapdfs_coordtransform[timpact[kk]]._kick_interpolatedThetasTrack
      
            
        #at each timpact compute apar_max
        apar_max=sdf_smooth.length(tdisrupt=sdf_pepper._tdisrupt-timpact[kk])*sdf_pepper._length_factor
        #print (apar_max)
        #considering the stream until apar_max, store xyzvxvyvz 
        for ii in range(len(apar_all)):
            if apar_all[ii] <= apar_max :
                apar.append(apar_all[ii])
                x.append(a[:,0][ii])
                y.append(a[:,1][ii])
                z.append(a[:,2][ii])
                vx.append(a[:,3][ii])
                vy.append(a[:,4][ii])
                vz.append(a[:,5][ii])
            
        x_full.append(np.array(x))
        y_full.append(np.array(y))
        z_full.append(np.array(z))
        vx_full.append(np.array(vx))
        vy_full.append(np.array(vy))
        vz_full.append(np.array(vz))
        apar_full.append(np.array(apar)*sdf_pepper._sigMeanSign) # _sigMeanSign = -/+ = trail/lead
        
    return (apar_full,x_full,y_full,z_full,vx_full,vy_full,vz_full)
    

def aparxv_stream_from_pkl(sampling=256,nchunks=16):
    
    '''
    compute apar,x,v from one or multiple pickle files
    '''
    
    apar=[]
    x_stream=[]
    y_stream=[]
    z_stream=[]
    vx_stream=[]
    vy_stream=[]
    vz_stream=[]
    timpact=[]
    
    sdf_smooth= pal5_util.setup_pal5model()
    
    if nchunks > 1 :
        
        for i in range(nchunks):
            with open('pkl_files/pal5pepper_{}sampling_MW2014_{}.pkl'.format(sampling,i),'rb') as savefile:
                    
                    print (sampling,i)

                    sdf_pepper= pickle.load(savefile,encoding='latin1')
                    ap,x,y,z,vx,vy,vz= aparxv_stream(sdf_smooth,sdf_pepper)
                    apar.extend(ap)
                    x_stream.extend(x)
                    y_stream.extend(y)
                    z_stream.extend(z)
                    vx_stream.extend(vx)
                    vy_stream.extend(vy)
                    vz_stream.extend(vz)
                    timpact.extend(sdf_pepper._timpact)
                    
    else :
        
        with open('pkl_files/pal5pepper_{}sampling_MW2014.pkl'.format(sampling),'rb') as savefile:

                    sdf_pepper= pickle.load(savefile,encoding='latin1')
                    ap,x,y,z,vx,vy,vz= aparxv_stream(sdf_smooth,sdf_pepper)
                    apar.extend(ap)
                    x_stream.extend(x)
                    y_stream.extend(y)
                    z_stream.extend(z)
                    vx_stream.extend(vx)
                    vy_stream.extend(vy)
                    vz_stream.extend(vz)
                    timpact.extend(sdf_pepper._timpact)
        
                
    return (timpact,apar,x_stream,y_stream,z_stream,vx_stream,vy_stream,vz_stream)
    
    
def compute_impact_parameters(timp,a,xs,ys,zs,nchunks=16,sampling_low=128,imp_fac=5.,Mmin=10**6.,rand_rotate=False):
    
    '''
    timp : timpacts
    a,xs,ys,zs : list of array, each array decribes the stream at that time, no of arrays = timpacts
    sampling_low : low timpact object on to which the impacts from high timpact case will be set
    imp_fac: X where bmax= X.r_s
    Mmin min mass above which all GMCs will be considered for impact
    rand_rotate : give the GMCs an ol' shaka shaka along phi
    
    '''
       
    #load the GMCs
    M,rs,coord=add_MCs(Mmin=Mmin,rand_rotate=rand_rotate)

    #integrate their orbits 5 Gyr back,
    t_age= np.linspace(0.,5.,1001)/bovy_conversion.time_in_Gyr(vo,ro)

    orbits=[]

    N=len(M)

    for ii in range(N):
    
        orbits.append(Orbit(coord[ii]).flip()) # flip flips the velocities for backwards integration
        orbits[ii].integrate(t_age,MWPotential2014)
        
    min_sep_matrix=np.empty([N,len(timp)])
    apar_matrix=np.empty([N,len(timp)])

    #compute min_sep of each MC
    for kk in range(len(timp)):
        for jj in range(N) :
            x_mc=orbits[jj].x(timp[kk])
            y_mc=orbits[jj].y(timp[kk])
            z_mc=orbits[jj].z(timp[kk])

            min_sep,apar_min=compute_min_separation(x_mc,y_mc,z_mc,a[kk],xs[kk],ys[kk],zs[kk])

            min_sep_matrix[jj,kk]=min_sep
            apar_matrix[jj,kk]=apar_min
            
            
    impactb=[]
    impact_angle=[]
    vx_mc=[]
    vy_mc=[]
    vz_mc=[]
    tmin=[]
    rs_mc=[]
    M_mc=[]
    impactMC_ind=[]

    if nchunks > 1 :
        with open('pkl_files/pal5pepper_{}sampling_MW2014.pkl'.format(sampling_low),'rb') as savefile:
            sdf_pepper_low= pickle.load(savefile,encoding='latin1')
                
        timpact_low=sdf_pepper_low._timpact
        
        c=0
        for ii in range(len(orbits)):

            bmax=imp_fac*rs[ii]/ro

            if min(min_sep_matrix[ii]) <= bmax :
                c+=1

                min_timpact_ind=np.argmin(min_sep_matrix[ii])

                impactMC_ind.append(ii)

                t_high=timp[min_timpact_ind]

                #round t_high to the nearest timpact in the low timpact sampling
                t_low=timpact_low[np.argmin(np.abs(timpact_low-t_high))]
                tmin.append(t_low)

                impactb.append(min_sep_matrix[ii,min_timpact_ind])
                impact_angle.append(apar_matrix[ii,min_timpact_ind]) # _sigMeanSign = -/+ = trail/lead

                rs_mc.append(rs[ii]/ro)
                M_mc.append(M[ii]/bovy_conversion.mass_in_msol(vo,ro))
                #flip velocities
                vx_mc.append(-orbits[ii].vx(t_high))
                vy_mc.append(-orbits[ii].vy(t_high))
                vz_mc.append(-orbits[ii].vz(t_high))

        #combine vx,vy,vz to v
        v_mc=np.c_[vx_mc,vy_mc,vz_mc]
        print ("The stream had %i impacts"%c)
        
    else :
        
        c=0
        for ii in range(len(orbits)):

            bmax=imp_fac*rs[ii]/ro

            if min(min_sep_matrix[ii]) <= bmax :
                c+=1

                min_timpact_ind=np.argmin(min_sep_matrix[ii])

                impactMC_ind.append(ii)

                t_imp_min=timp[min_timpact_ind]
                tmin.append(t_imp_min)

                impactb.append(min_sep_matrix[ii,min_timpact_ind])
                impact_angle.append(apar_matrix[ii,min_timpact_ind]) # _sigMeanSign = -/+ = trail/lead

                rs_mc.append(rs[ii]/ro)
                M_mc.append(M[ii]/bovy_conversion.mass_in_msol(vo,ro))
                #flip velocities
                vx_mc.append(-orbits[ii].vx(t_imp_min))
                vy_mc.append(-orbits[ii].vy(t_imp_min))
                vz_mc.append(-orbits[ii].vz(t_imp_min))

        #combine vx,vy,vz to v
        v_mc=np.c_[vx_mc,vy_mc,vz_mc]
        print ("The stream had %i impacts"%c)
        
        
    
    return (impactMC_ind,M_mc,rs_mc,v_mc,impactb,impact_angle,tmin)