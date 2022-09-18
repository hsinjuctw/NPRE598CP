import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import Library.bearth as bearth
import Library.coord as coord
import Library.particle as particle

def main():
    ##################################
    ##### EARTH'S MAGNETIC FIELD #####
    ##################################
    RE = 6.371e6 # Earth radius [m]
    Nx = 21
    Ny = 21
    Nz = 21
    #r     = np.linspace( RE,   20.*RE, Nr )
    #theta = np.linspace( 0.,    np.pi, Nt )
    #phi   = np.linspace( 0., 2.*np.pi, Np )
    X = np.linspace( -20.*RE, 20.*RE, Nx)
    Y = np.linspace( -20.*RE, 20.*RE, Ny)
    Z = np.linspace( -20.*RE, 20.*RE, Nz)
    
    Bx,By,Bz,Bmag = bearth.getEarthDipole(X,Y,Z)
    
    refine  = 3.
    Xrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Nx-1)+1))
    Yrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Ny-1)+1))
    Zrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Nz-1)+1))
    # interpolated b field array initialization
    Bxi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
    Byi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
    Bzi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
    for i in range(0,int(refine*(Nx-1)+1)):
        for j in range(0,int(refine*(Ny-1)+1)):
            for k in range(0,int(refine*(Nz-1)+1)):
                if abs(Xrefine[i])<RE and abs(Yrefine[j])<RE and abs(Zrefine[k])<RE:
                    Bxi[i,j,k] = np.nan
                    Byi[i,j,k] = np.nan
                    Bzi[i,j,k] = np.nan
                else:
                    xp = Xrefine[i]
                    yp = Yrefine[j]
                    zp = Zrefine[k]
                    Bxi[i,j,k],Byi[i,j,k],Bzi[i,j,k] = particle.interp3D(X,Y,Z,Bx,By,Bz,xp,yp,zp)
    
    plt.figure(1)
    XX,ZZ = np.meshgrid(X,Z)
    plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bx[:,10,:])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('Bx in XZ plane [T] - Earth magnetic field')

    plt.figure(2)
    XXrefine,ZZrefine = np.meshgrid(Xrefine,Zrefine)
    plt.contourf(np.transpose(XXrefine)/RE,np.transpose(ZZrefine)/RE, Bxi[:,int(10*refine),:])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('Interpolated Bx in XZ plane [T] - Earth magnetic field')
    plt.show()

if __name__ == '__main__':
   main()