import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import Library.bearth as bearth
import Library.coord as coord
from matplotlib import ticker#, cm


def main():
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
    
    #Bx,By,Bz,Bmag = bearth.getEarthDipole(X,Y,Z)
    Bx,By,Bz,Bmag = bearth.getEarthDipoleCSC(X,Y,Z)

    # print(np.nanmax(Bmag))
    # print(np.nanmin(Bmag[:,:,25]))

    # X0  = 10.*RE
    # Y0  = 0.
    # Z0  = 0.
    # X0  = .35*RE
    # Y0  = 0.
    # Z0  = -.94*RE
    Lshell = 10.
    MLT    = 12.
    dX  = .01*RE
    # dY  = .01*RE
    dZ  = .01*RE
    # rdtheta = np.pi/1800.

    FLx0,FLy0,FLz0 = bearth.dipoleFieldline2D(5.,12.,dX,dZ)
    FLx1,FLy1,FLz1 = bearth.dipoleFieldline2D(10.,12.,dX,dZ)
    FLx2,FLy2,FLz2 = bearth.dipoleFieldline2D(15.,12.,dX,dZ)
    FLx3,FLy3,FLz3 = bearth.dipoleFieldline2D(5.,0.,dX,dZ)
    FLx4,FLy4,FLz4 = bearth.dipoleFieldline2D(10.,0.,dX,dZ)
    FLx5,FLy5,FLz5 = bearth.dipoleFieldline2D(15.,0.,dX,dZ)

    cmap = 'PRGn'
    # cmap = 'RdBu'
    # plt.figure(1)
    fig, ax1 = plt.subplots()
    XX,ZZ = np.meshgrid(X,Z)
    plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bmag[:,50,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
    ax1.add_artist(Circle((0,0), 1, color='b'))
    plt.plot( FLx0/RE,FLz0/RE, 'k-' )
    plt.plot( FLx1/RE,FLz1/RE, 'k-' )
    plt.plot( FLx2/RE,FLz2/RE, 'k-' )
    plt.plot( FLx3/RE,FLz3/RE, 'k-' )
    plt.plot( FLx4/RE,FLz4/RE, 'k-' )
    plt.plot( FLx5/RE,FLz5/RE, 'k-' )
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in XZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XZ.png',dpi=150)
    plt.show()
    
    # plt.figure(2)
    fig, ax2 = plt.subplots()
    XX,YY = np.meshgrid(X,Y)
    plt.contourf(np.transpose(XX)/RE,np.transpose(YY)/RE, Bmag[:,:,50],locator=ticker.LogLocator(),cmap=cmap)
    ax2.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Y$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in XY plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XY.png',dpi=150)
    plt.show()
    
    # plt.figure(3)
    fig, ax3 = plt.subplots()
    YY,ZZ = np.meshgrid(Y,Z)
    plt.contourf(np.transpose(YY)/RE,np.transpose(ZZ)/RE, Bmag[50,:,:],locator=ticker.LogLocator(),cmap=cmap)
    ax3.add_artist(Circle((0,0), 1, color='b'))
    plt.xlabel('$Y$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in YZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_YZ.png',dpi=150)
    plt.show()

if __name__ == '__main__':
   main()