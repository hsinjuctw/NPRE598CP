import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import Library.bearth as bearth
import Library.coord as coord
from matplotlib import ticker#, cm
import time

def main():
    # record start time
    start = time.time()
    ##################################
    ##### EARTH'S MAGNETIC FIELD #####
    ##################################
    RE = 6.371e6 # Earth radius [m]
    Nx = 101
    Ny = 101
    Nz = 101
    #r     = np.linspace( RE,   20.*RE, Nr )
    #theta = np.linspace( 0.,    np.pi, Nt )
    #phi   = np.linspace( 0., 2.*np.pi, Np )
    X = np.linspace( -20.*RE, 20.*RE, Nx)
    Y = np.linspace( -20.*RE, 20.*RE, Ny)
    Z = np.linspace( -20.*RE, 20.*RE, Nz)
    
    Bx,By,Bz,Bmag = bearth.getEarthDipole(X,Y,Z)
    #Bx,By,Bz,Bmag = bearth.getEarthDipoleCSC(X,Y,Z)


    Lshell = 10.
    MLT    = 12.
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

    # record end time
    end = time.time()
    print('The time of execution of above program is :',(end-start) * 10**3, ' ms')

    #################
    ##### PLOTS #####
    #################
    cmap = 'PRGn'

    plt.rc('legend',fontsize=12)
    plt.rc('axes',labelsize=12)
    plt.rc('xtick',labelsize=12)
    plt.rc('ytick',labelsize=12)
    # fig, ax1 = plt.subplots(figsize=(4.8, 4.8))
    fig, ax1 = plt.subplots()
    XX,ZZ = np.meshgrid(X,Z)
    plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bmag[:,50,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    ax1.add_artist(Circle((0,0), 1, color='b'))
    plt.plot( FLx0/RE,FLz0/RE, 'b-', linewidth=.5 )
    plt.plot( FLx1/RE,FLz1/RE, 'b-', linewidth=.5 )
    plt.plot( FLx2/RE,FLz2/RE, 'b-', linewidth=.5 )
    plt.plot( FLx3/RE,FLz3/RE, 'b-', linewidth=.5 )
    plt.plot( FLx4/RE,FLz4/RE, 'b-', linewidth=.5 )
    plt.plot( FLx5/RE,FLz5/RE, 'b-', linewidth=.5 )
    plt.plot( FLx6/RE,FLz6/RE, 'b-', linewidth=.5 )
    plt.plot( FLx7/RE,FLz7/RE, 'b-', linewidth=.5 )
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    cbar = plt.colorbar(format='%.0e')
    cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax1.set_aspect('equal', 'box')
    # plt.legend(loc=3)
    plt.tight_layout()
    # plt.title('B-field magnitude in XZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XZ.png',dpi=150)
    plt.show()

    fig, ax2 = plt.subplots(figsize=(4.8, 4.8))
    XX,YY = np.meshgrid(X,Y)
    plt.contourf(np.transpose(XX)/RE,np.transpose(YY)/RE, Bmag[:,:,50],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    ax2.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Y$ [$R_E$]')
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax2.set_aspect('equal', 'box')
    # plt.legend(loc=3)
    plt.tight_layout()
    # plt.title('B-field magnitude in XY plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XY.png',dpi=150)
    plt.show()
    
    # plt.figure(3)
    fig, ax3 = plt.subplots()
    YY,ZZ = np.meshgrid(Y,Z)
    plt.contourf(np.transpose(YY)/RE,np.transpose(ZZ)/RE, Bmag[50,:,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    ax3.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$Y$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    # cbar = plt.colorbar(format='%.0e')
    # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
    ax3.set_aspect('equal', 'box')
    # plt.legend(loc=3)
    plt.tight_layout()
    # plt.title('B-field magnitude in YZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_YZ.png',dpi=150)
    plt.show()

if __name__ == '__main__':
   main()