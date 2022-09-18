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
    x = np.linspace( -20.*RE, 20.*RE, Nx)
    y = np.linspace( -20.*RE, 20.*RE, Ny)
    z = np.linspace( -20.*RE, 20.*RE, Nz)
    
    Bmagc = np.empty((Nx,Ny,Nz))
    Bxc   = np.empty((Nx,Ny,Nz))
    Byc   = np.empty((Nx,Ny,Nz))
    Bzc   = np.empty((Nx,Ny,Nz))
    
    for i in range(0,Nx):
        for j in range(0,Ny):
            for k in range(0,Nz):
                #if abs(x[i])<RE and abs(y[j])<RE and abs(z[k])<RE:
                if abs(x[i]/RE)==0. and abs(y[j]/RE)==0. and abs(z[k]/RE)==0.:
                    Bmagc[i,j,k] = np.nan
                    Bxc[i,j,k]   = np.nan
                    Byc[i,j,k]   = np.nan
                    Bzc[i,j,k]   = np.nan
                else:
                    Bxc[i,j,k],Byc[i,j,k],Bzc[i,j,k] = bearth.dipoleEarth(x[i],y[j],z[k])
                    Bmagc[i,j,k] = np.sqrt( Bxc[i,j,k]**2 + Byc[i,j,k]**2 + Bzc[i,j,k]**2 )
    #xp = 2.5*RE
    #yp = 2.2*RE
    #zp = 3.2*RE
    #print(particle.interp3Dalt(x,y,z,Bxc,Byc,Bzc,xp,yp,zp))
    #print(bearth.dipoleEarth(xp,yp,zp))
    
    refine  = 2.
    xrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Nx-1)+1))
    yrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Ny-1)+1))
    zrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Nz-1)+1))
    # interpolated b field array initialization
    Bxi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
    Byi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
    Bzi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
    #j     = 50
    for i in range(0,int(refine*(Nx-1)+1)):
        for j in range(0,int(refine*(Ny-1)+1)):
            # xp = xrefine[i]
            # zp = zrefine[k]
            # Bxi[i,j,k],Byi[i,j,k],Bzi[i,j,k] = particle.interp2D(x,z,Bxc[:,10,:],Byc[:,10,:],Bzc[:,10,:],xp,zp)
            for k in range(0,int(refine*(Nz-1)+1)):
                if abs(xrefine[i])<RE and abs(yrefine[j])<RE and abs(zrefine[k])<RE:
                    Bxi[i,j,k] = np.nan
                    Byi[i,j,k] = np.nan
                    Bzi[i,j,k] = np.nan
                else:
                    if i%refine == 0 and j%refine == 0 and k%refine == 0:
                        Bxi[i,j,k] = Bxc[int(i/refine),int(j/refine),int(k/refine)]
                        Byi[i,j,k] = Byc[int(i/refine),int(j/refine),int(k/refine)]
                        Bzi[i,j,k] = Bzc[int(i/refine),int(j/refine),int(k/refine)]
                        #print('hi')
                    else:
                        xp = xrefine[i]
                        yp = yrefine[j]
                        zp = zrefine[k]
                        Bxi[i,j,k],Byi[i,j,k],Bzi[i,j,k] = particle.interp3D(x,y,z,Bxc,Byc,Bzc,xp,yp,zp)
    
    plt.figure(1)
    xx,zz = np.meshgrid(x,z)
    plt.contourf(np.transpose(xx)/RE,np.transpose(zz)/RE, Bxc[:,10,:])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('Bx in XZ plane [T] - Earth magnetic field')

    plt.figure(2)
    xxrefine,zzrefine = np.meshgrid(xrefine,zrefine)
    plt.contourf(np.transpose(xxrefine)/RE,np.transpose(zzrefine)/RE, Bxi[:,int(10*refine),:])
    plt.xlabel('$X$ [$R_E$]')
    plt.ylabel('$Z$ [$R_E$]')
    plt.axis('equal')
    plt.colorbar()
    plt.title('Interpolated Bx in XZ plane [T] - Earth magnetic field')
    plt.show()

if __name__ == '__main__':
   main()