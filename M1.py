import numpy as np
import time as t
import matplotlib.pyplot as plt
import Library.bearth as bearth
import Library.particle as particle
import Library.coord as coord
from matplotlib.patches import Circle
from matplotlib import ticker

qe = 1.60217662e-19
me = 9.10938356e-31
mp = 1.6726219e-27
RE = 6.371e6 # Earth radius [m]

def main():
    # record start time
    start = t.time()

    ##################################
    ##### EARTH'S MAGNETIC FIELD #####
    ##################################
    Nx = 101
    Ny = 101
    Nz = 101

    X = np.linspace( -20.*RE, 20.*RE, Nx)
    Y = np.linspace( -20.*RE, 20.*RE, Ny)
    Z = np.linspace( -20.*RE, 20.*RE, Nz)
    
    Bx,By,Bz,Bmag = bearth.getEarthDipole(X,Y,Z)

    ###############################
    ##### PARTICLE TRAJECTORY #####
    ###############################
    fcorrection = True # flag for frequency correction

    # charge [C], mass [kg]
    Q = -qe
    M = me

    # # cyclotron frequency [rad/s]
    # omega_c = np.abs(Q)*Bmag/M
    # # cyclotron period [s]
    # Tc = 2.*np.pi/omega_c
    # # magnitude of v_perp
    # v0 = np.sqrt(vx0**2+vy0**2)
    # # Larmor radius
    # r_L = v0/omega_c

    Tev = 1.e4             # particle temperature in eV
    T = particle.eV2K(Tev) # particle temperature in K
    v = particle.E2v(Tev,M)
    Np = 10000            # particle count

    Lr = 19.*RE
    Lp = 2.*np.pi
    Rr = np.random.rand(Np)*Lr+RE
    Rp = np.random.rand(Np)*Lp
    # particle [i] is at (Rr[i],np.pi/2.,Rp[i])
    Loc = np.empty((Np,3)) # particle location in Cartesian

    Vx = np.linspace(-3.*v,3.*v,100) # arbitrary end points
    Xp = np.random.rand(Np)
    Yp = np.random.rand(Np)
    Zp = np.random.rand(Np)

    Vpx = particle.invMaxwellianCDF1D(M,T,Vx,Xp)
    Vpy = particle.invMaxwellianCDF1D(M,T,Vx,Yp)
    Vpz = particle.invMaxwellianCDF1D(M,T,Vx,Zp)

    # time grids [s]
    timetot = 10
    N_per_sec = 100
    time = np.linspace(0.,timetot,N_per_sec*timetot+1)
    dt = time[1]-time[0]

    BB = np.empty((Np,len(time),6))

    for n in range(0,Np):
        Loc[n,0],Loc[n,1],Loc[n,2] = coord.sph2car(Rr[n],np.pi/2.,Rp[n])
        x0 = Loc[n,0]
        y0 = Loc[n,1]
        z0 = Loc[n,2]
        # B field [nT]
        bx, by, bz = bearth.dipoleEarth(x0,y0,z0)
        b = np.sqrt(bx**2+by**2+bz**2)
        # initial velocity [m/s]
        vx0 = Vpx[n]
        vy0 = Vpy[n]
        vz0 = Vpz[n]
        # initial state vector
        BB0 = np.array((x0,y0,z0,vx0,vy0,vz0))
        # parameters
        params = np.array((dt, Q/M*dt/2.))
        # Boris integration
        BB[n,:,:] = particle.dirBorisBunemann(time, BB0, params, fcorrection)
    
    # record end time
    end = t.time()
    print('The time spent for solving the particles is ',(end-start) * 10**3, ' ms')

    ##############################
    ##### FIELD LINE TRACING #####
    ##############################
    # Lshell = 10.
    # MLT    = 12.
    dX  = .01*RE
    dZ  = .01*RE

    FLx0,FLy0,FLz0 = bearth.dipoleFieldline2D(5.,12.,dX,dZ)
    FLx1,FLy1,FLz1 = bearth.dipoleFieldline2D(5.,0.,dX,dZ)
    FLx2,FLy2,FLz2 = bearth.dipoleFieldline2D(10.,12.,dX,dZ)
    FLx3,FLy3,FLz3 = bearth.dipoleFieldline2D(10.,0.,dX,dZ)
    FLx4,FLy4,FLz4 = bearth.dipoleFieldline2D(15.,12.,dX,dZ)
    FLx5,FLy5,FLz5 = bearth.dipoleFieldline2D(15.,0.,dX,dZ)
    FLx6,FLy6,FLz6 = bearth.dipoleFieldline2D(20.,12.,dX,dZ)
    FLx7,FLy7,FLz7 = bearth.dipoleFieldline2D(20.,0.,dX,dZ)

    #################
    ##### PLOTS #####
    #################
    cmap = 'PRGn'

    plt.rc('legend',fontsize=12)
    plt.rc('axes',labelsize=12)
    plt.rc('xtick',labelsize=12)
    plt.rc('ytick',labelsize=12)

    fig = plt.figure(figsize=(8, 7.6))
    ax = fig.add_subplot(111)
    plt.plot( np.transpose(BB[:,:,0])/RE, np.transpose(BB[:,:,2])/RE )
    XX,ZZ = np.meshgrid(X,Z)
    plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bmag[:,50,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    plt.plot( FLx0/RE,FLz0/RE, 'b-', linewidth=.5 )
    plt.plot( FLx1/RE,FLz1/RE, 'b-', linewidth=.5 )
    plt.plot( FLx2/RE,FLz2/RE, 'b-', linewidth=.5 )
    plt.plot( FLx3/RE,FLz3/RE, 'b-', linewidth=.5 )
    plt.plot( FLx4/RE,FLz4/RE, 'b-', linewidth=.5 )
    plt.plot( FLx5/RE,FLz5/RE, 'b-', linewidth=.5 )
    plt.plot( FLx6/RE,FLz6/RE, 'b-', linewidth=.5 )
    plt.plot( FLx7/RE,FLz7/RE, 'b-', linewidth=.5 )
    ax.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$x$ [$R_E$]')
    plt.ylabel('$z$ [$R_E$]')
    plt.xlim([-20.,20.])
    plt.ylim([-20.,20.])
    # plt.title('Trajectory projection on $XZ$-plane')
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax.set_aspect('equal', 'box')
    plt.savefig('M1_XZ.png')

    fig = plt.figure(figsize=(8, 7.6))
    ax = fig.add_subplot(111)
    plt.plot( np.transpose(BB[:,:,0])/RE, np.transpose(BB[:,:,1])/RE )
    XX,YY = np.meshgrid(X,Y)
    plt.contourf(np.transpose(XX)/RE,np.transpose(YY)/RE, Bmag[:,:,50],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    ax.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$x$ [$R_E$]')
    plt.ylabel('$y$ [$R_E$]')
    plt.xlim([-20.,20.])
    plt.ylim([-20.,20.])
    # plt.title('Trajectory projection on $XY$-plane')
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax.set_aspect('equal', 'box')
    plt.savefig('M1_XY.png')

    fig = plt.figure(figsize=(8, 7.6))
    ax = fig.add_subplot(111)
    plt.plot( np.transpose(BB[:,:,1])/RE, np.transpose(BB[:,:,2])/RE )
    YY,ZZ = np.meshgrid(Y,Z)
    plt.contourf(np.transpose(YY)/RE,np.transpose(ZZ)/RE, Bmag[50,:,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    ax.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$y$ [$R_E$]')
    plt.ylabel('$z$ [$R_E$]')
    plt.xlim([-20.,20.])
    plt.ylim([-20.,20.])
    # plt.title('Trajectory projection on $YZ$-plane')
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax.set_aspect('equal', 'box')
    plt.savefig('M1_YZ.png')
    plt.show()

if __name__ == '__main__':
   main()