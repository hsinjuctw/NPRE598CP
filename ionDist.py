import numpy as np
import matplotlib.pyplot as plt
import Library.bearth as bearth
import Library.particle as particle
import Library.coord as coord
from matplotlib.patches import Circle

qe = 1.60217662e-19
me = 9.10938356e-31
mp = 1.6726219e-27

def main():
    RE = 6.371e6 # Earth radius [m]

    ##################################
    ##### EARTH'S MAGNETIC FIELD #####
    ##################################
    # SPHERICAL CASE
    # Nr = 101
    # Nt = 101

    # R     = np.linspace( RE,   20.*RE, Nr )
    # Theta = np.linspace( 0., 2.*np.pi, Nt )

    # Br,Btheta,Bphi = bearth.getEarthDipoleSph(R,Theta)
    # Bmag = Br**2+Btheta**2+Bphi**2
    # print(np.nanmax(Bmag))
    # print(np.nanmin(Bmag))

    ###############################
    ##### PARTICLE TRAJECTORY #####
    ###############################
    fcorrection = True # flag for frequency correction

    # charge [C], mass [kg]
    Q = -qe
    M = me

    # # # cyclotron frequency [rad/s]
    # # omega_c = np.abs(Q)*Bmag/M
    # # # cyclotron period [s]
    # # Tc = 2.*np.pi/omega_c
    # # # magnitude of v_perp
    # # v0 = np.sqrt(vx0**2+vy0**2)
    # # # Larmor radius
    # # r_L = v0/omega_c

    Tev = 1.e4             # particle temperature in eV
    T = particle.eV2K(Tev) # particle temperature in K
    v = particle.E2v(Tev,M)
    Np = 10000            # particle count

    Lr = 19.*RE
    Lp = 2.*np.pi
    # Rr = np.random.rand(Np)*Lr+RE
    Rr = np.ones(Np)*10.*RE
    Rp = np.random.rand(Np)*Lp
    # particle [i] is at (Rr[i],np.pi/2.,Rp[i])
    Loc = np.empty((Np,3)) # particle location in Cartesian

    Vx = np.linspace(-3.*v,3.*v,100) # arbitrary end points
    print(particle.maxwellianCDF1D(M,3.*v,T))
    noP = np.linspace(0,Np,Np)
    Xp = np.random.rand(Np)
    Yp = np.random.rand(Np)
    Zp = np.random.rand(Np)

    Vpx = particle.invMaxwellianCDF1D(M,T,Vx,Xp)
    Vpy = particle.invMaxwellianCDF1D(M,T,Vx,Yp)
    Vpz = particle.invMaxwellianCDF1D(M,T,Vx,Zp)

    ###################################
    ##### PLOT SPEED DISTRUBUTION #####
    ###################################
    plt.rc('legend',fontsize=12)
    plt.rc('axes',labelsize=12)
    plt.rc('xtick',labelsize=12)
    plt.rc('ytick',labelsize=12)
    plt.figure(figsize=(6, 3))
    plt.hist(np.sqrt(Vpx**2+Vpy**2+Vpz**2),color = 'green',alpha=.5,bins=100)
    plt.xlabel('speed [m/s]')
    plt.ylabel('count')
    plt.tight_layout()
    plt.savefig('mb_reconstruct.png',dpi=150)
    plt.show()
    # fcorrection = True

    # # time grids [s]
    # timetot = 10
    # N_per_sec = 100
    # time = np.linspace(0.,timetot,N_per_sec*timetot+1)
    # dt = time[1]-time[0]

    # BB = np.empty((Np,len(time),6))
    # # BB = np.empty((Np,len(time),10))
    # # mu = np.empty((Np,len(time)))
    # E = np.empty((Np,len(time)))

    # for n in range(0,Np):
    #     Loc[n,0],Loc[n,1],Loc[n,2] = coord.sph2car(Rr[n],np.pi/2.,Rp[n])
    #     x0 = Loc[n,0]
    #     y0 = Loc[n,1]
    #     z0 = Loc[n,2]
    #     #x0,y0,z0 = coord.sph2car(2.5*RE,np.pi/8.,0.)
    #     # B field [nT]
    #     bx, by, bz = bearth.dipoleEarth(x0,y0,z0)
    #     b = np.sqrt(bx**2+by**2+bz**2)
    #     # initial velocity [m/s]
    #     vx0 = Vpx[n]
    #     vy0 = Vpy[n]
    #     vz0 = Vpz[n]
    #     # initial state vector
    #     BB0 = np.array((x0,y0,z0,vx0,vy0,vz0))
    #     # parameters
    #     params = np.array((dt, Q/M*dt/2.))
    #     # Boris integration
    #     BB[n,:,:] = particle.dirBorisBunemann(time, BB0, params, fcorrection)
    #     for t in range(0,len(time)):
    #         E[n,t] = .5*M*(BB[n,t,3]**2+BB[n,t,4]**2+BB[n,t,5]**2)
    # #     #     mu[n,t] = .5*particle.lorentz(np.sqrt(BB[n,t,3]**2+BB[n,t,4]**2+BB[n,t,5]**2))*M*(BB[n,t,7]**2+BB[n,t,8]**2+BB[n,t,9]**2)/BB[n,t,6]

    # fig = plt.figure(5)
    # ax = fig.add_subplot(111)
    # plt.plot( time, np.transpose(E))
    # plt.xlabel('time [s]')
    # # plt.ylabel('$\mu$')
    # plt.ylabel('$E_k$')
    # # plt.title('Magnetic Moment')
    # plt.title('Kinetic Energy')
    # plt.legend(loc=3)

    # fig = plt.figure(0)
    # ax = fig.add_subplot(111)
    # plt.plot( noP, Vpx, label='Initial Vx' )
    # plt.plot( noP, Vpy, label='Initial Vy' )
    # plt.plot( noP, Vpz, label='Initial Vz' )
    # plt.xlabel('particle')
    # plt.ylabel('$v$ [m/s]')
    # plt.title('Initial Velocity')
    # plt.legend(loc=3)

    # fig = plt.figure(1)
    # ax = fig.add_subplot(111)
    # plt.plot( np.transpose(BB[:,:,0])/RE, np.transpose(BB[:,:,2])/RE )
    # ax.add_artist(Circle((0,0), 1, color='b'))
    # plt.xlabel('$x$ [$R_E$]')
    # plt.ylabel('$z$ [$R_E$]')
    # plt.title('Trajectory projection on $XZ$-plane')
    # plt.axis('equal')
    # plt.legend(loc=3)
    # plt.savefig('boris_XZ.png')

    # fig = plt.figure(2)
    # ax = fig.add_subplot(111)
    # plt.plot( np.transpose(BB[:,:,0])/RE, np.transpose(BB[:,:,1])/RE )
    # ax.add_artist(Circle((0,0), 1, color='b'))
    # plt.xlabel('$x$ [$R_E$]')
    # plt.ylabel('$y$ [$R_E$]')
    # plt.axis('equal')
    # plt.legend(loc=3)
    # plt.savefig('boris_XY.png')
    # fig = plt.figure(3)
    # ax = fig.add_subplot(111)
    # plt.plot( np.transpose(BB[:,:,1])/RE, np.transpose(BB[:,:,2])/RE )
    # ax.add_artist(Circle((0,0), 1, color='b'))
    # plt.xlabel('$y$ [$R_E$]')
    # plt.ylabel('$z$ [$R_E$]')
    # plt.axis('equal')
    # plt.legend(loc=3)
    # plt.savefig('boris_YZ.png')
    # # plt.figure(4)
    # # plt.plot( time, X[:,0]/RE, 'r-', label='$x(t)$')
    # # plt.plot( time, X[:,1]/RE, 'g-', label='$y(t)$')
    # # plt.plot( time, X[:,2]/RE, 'b-', label='$z(t)$')
    # # plt.xlabel('$t$ [s]')
    # # plt.ylabel('location [$R_E$]')
    # # plt.legend(loc=3)
    # # plt.savefig('boris_time.png')
    # plt.show()

if __name__ == '__main__':
   main()