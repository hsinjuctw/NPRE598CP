import numpy as np
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
    x0 = -10.*RE
    y0 = 0.
    z0 = 0.
    #x0,y0,z0 = coord.sph2car(2.5*RE,np.pi/8.,0.)
    # B field [nT]
    bx, by, bz = bearth.dipoleEarth(x0,y0,z0)
    b = np.sqrt(bx**2+by**2+bz**2)
    # cyclotron frequency [rad/s]
    omega_c = np.abs(Q)*b/M
    # cyclotron period [s]
    # Tc = 2.*np.pi/omega_c
    # initial velocity [m/s]
    vx0 = -1.e7
    vy0 =  1.e7
    vz0 =  1.e7
    # magnitude of v_perp
    # v0 = np.sqrt(vx0**2+vy0**2)
    # Larmor radius
    # r_L = v0/omega_c
    # initial position [m]
    #x0 += r_L
    # initial state vector
    X0 = np.array((x0,y0,z0,vx0,vy0,vz0))
    # time grids [s]
    timetot = 250
    N_per_sec = 100
    time = np.linspace(0.,timetot,N_per_sec*timetot+1)
    dt = time[1]-time[0]
    # parameters
    params = np.array((dt, Q/M*dt/2.))
    # Boris integration
    BB = particle.dirBorisBunemann(time, X0, params, fcorrection)
    
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

    # fig = plt.figure(1)

    fig = plt.figure(figsize=(4.8, 4.8))
    plt.rc('legend',fontsize=12)
    plt.rc('axes',labelsize=12)
    plt.rc('xtick',labelsize=12)
    plt.rc('ytick',labelsize=12)

    ax = fig.add_subplot(111)
    XX,ZZ = np.meshgrid(X,Z)
    plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bmag[:,50,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    plt.plot( BB[:,0]/RE, BB[:,2]/RE, 'k-', label='Trajectory projection on $xz$-plane', linewidth=1. )
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
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax.set_aspect('equal', 'box')
    plt.legend(loc=3)
    plt.tight_layout()
    plt.savefig('boris_XZ.png')


    # fig = plt.figure(2)
    fig = plt.figure(figsize=(4.8, 4.8))
    ax = fig.add_subplot(111)
    XX,YY = np.meshgrid(X,Y)
    plt.contourf(np.transpose(XX)/RE,np.transpose(YY)/RE, Bmag[:,:,50],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    plt.plot( BB[:,0]/RE, BB[:,1]/RE, 'k-', label='Trajectory projection on $xy$-plane', linewidth=1.  )
    ax.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$x$ [$R_E$]')
    plt.ylabel('$y$ [$R_E$]')
    plt.xlim([-20,20])
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax.set_aspect('equal', 'box')
    plt.legend(loc=3)
    plt.tight_layout()
    plt.savefig('boris_XY.png')

    # fig = plt.figure(3)
    # ax = fig.add_subplot(111)
    # plt.plot( BB[:,1]/RE, BB[:,2]/RE, 'k-', label='Trajectory projection on $yz$-plane' )
    # ax.add_artist(Circle((0,0), 1, color='b'))
    # plt.xlabel('$y$ [$R_E$]')
    # plt.ylabel('$z$ [$R_E$]')
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    # ax.set_aspect('equal', 'box')
    # plt.legend(loc=3)
    # plt.savefig('boris_YZ.png')

    plt.figure(4)
    plt.plot( time, BB[:,0]/RE, 'r-', label='$x(t)$')
    plt.plot( time, BB[:,1]/RE, 'g-', label='$y(t)$')
    plt.plot( time, BB[:,2]/RE, 'b-', label='$z(t)$')
    plt.xlabel('$t$ [s]')
    plt.ylabel('location [$R_E$]')
    plt.legend(loc=3)
    plt.savefig('boris_time.png')
    plt.show()

if __name__ == '__main__':
   main()