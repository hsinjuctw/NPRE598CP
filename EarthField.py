import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import Library.bearth as bearth
import Library.coord as coord


def main():
    ##################################
    ##### EARTH'S MAGNETIC FIELD #####
    ##################################
    RE = 6.371e6 # Earth radius [m]
    Nx = 51
    Ny = 51
    Nz = 51
    #r     = np.linspace( RE,   20.*RE, Nr )
    #theta = np.linspace( 0.,    np.pi, Nt )
    #phi   = np.linspace( 0., 2.*np.pi, Np )
    X = np.linspace( -20.*RE, 20.*RE, Nx)
    Y = np.linspace( -20.*RE, 20.*RE, Ny)
    Z = np.linspace( -20.*RE, 20.*RE, Nz)
    
    Bx,By,Bz,Bmag = bearth.getEarthDipole(X,Y,Z)

    plt.figure(1)
    XX,ZZ = np.meshgrid(X,Z)
    plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bmag[:,25,:])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in XZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XZ.png',dpi=150)
    plt.show()
    
    plt.figure(2)
    XX,YY = np.meshgrid(X,Y)
    plt.contourf(np.transpose(XX)/RE,np.transpose(YY)/RE, Bmag[:,:,25])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Y$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in XY plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XY.png',dpi=150)
    plt.show()
    
    plt.figure(3)
    YY,ZZ = np.meshgrid(Y,Z)
    plt.contourf(np.transpose(YY)/RE,np.transpose(ZZ)/RE, Bmag[25,:,:])
    plt.xlabel('$Y$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in YZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_YZ.png',dpi=150)
    plt.show()

if __name__ == '__main__':
   main()