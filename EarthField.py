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
    x = np.linspace( -20.*RE, 20.*RE, Nx)
    y = np.linspace( -20.*RE, 20.*RE, Ny)
    z = np.linspace( -20.*RE, 20.*RE, Nz)
    
    Bmag = np.empty((Nx,Ny,Nz))
    Bx   = np.empty((Nx,Ny,Nz))
    By   = np.empty((Nx,Ny,Nz))
    Bz   = np.empty((Nx,Ny,Nz))
    
    # Bmagc = np.empty((Nx,Ny,Nz))
    # Bxc   = np.empty((Nx,Ny,Nz))
    # Byc   = np.empty((Nx,Ny,Nz))
    # Bzc   = np.empty((Nx,Ny,Nz))
    
    # for i in range(0,Nx):
    #     for j in range(0,Ny):
    #         for k in range(0,Nz):
    #             if abs(x[i])<RE and abs(y[j])<RE and abs(z[k])<RE:
    #                 Bmagc[i,j,k] = np.nan
    #                 Bxc[i,j,k]   = np.nan
    #                 Byc[i,j,k]   = np.nan
    #                 Bzc[i,j,k]   = np.nan
    #             else:
    #                 Bxc[i,j,k],Byc[i,j,k],Bzc[i,j,k] = bearth.dipoleEarth(x[i],y[j],z[k])
    #                 Bmagc[i,j,k] = np.sqrt( Bxc[i,j,k]**2 + Byc[i,j,k]**2 + Bzc[i,j,k]**2 )
    
    # calculate in spherical
    for i in range(0,Nx):
        for j in range(0,Ny):
            for k in range(0,Nz):
                if abs(x[i])<RE and abs(y[j])<RE and abs(z[k])<RE:
                    Bmag[i,j,k] = np.nan
                    Bx[i,j,k]   = np.nan
                    By[i,j,k]   = np.nan
                    Bz[i,j,k]   = np.nan
                else:
                    r,theta,phi = coord.car2sph( x[i],y[j],z[k] )
                    Br,Bt,Bp    = bearth.dipoleEarthSph( r,theta,phi )
                    Bx[i,j,k],By[i,j,k],Bz[i,j,k] = coord.sph2carV(Br,Bt,Bp,r,theta,phi)
                    Bmag[i,j,k] = np.sqrt( Bx[i,j,k]**2 + By[i,j,k]**2 + Bz[i,j,k]**2 )
    
    plt.figure(1)
    xx,zz = np.meshgrid(x,z)
    plt.contourf(np.transpose(xx)/RE,np.transpose(zz)/RE, Bmag[:,25,:])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in XZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XZ.png',dpi=150)
    plt.show()
    
    plt.figure(2)
    xx,yy = np.meshgrid(x,y)
    plt.contourf(np.transpose(xx)/RE,np.transpose(yy)/RE, Bmag[:,:,25])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Y$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in XY plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_XY.png',dpi=150)
    plt.show()
    
    plt.figure(3)
    yy,zz = np.meshgrid(y,z)
    plt.contourf(np.transpose(yy)/RE,np.transpose(zz)/RE, Bmag[25,:,:])
    plt.xlabel('$Y$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('B-field magnitude in YZ plane [T] - Earth magnetic field')
    plt.savefig('dipoleEarth_YZ.png',dpi=150)
    plt.show()

if __name__ == '__main__':
   main()