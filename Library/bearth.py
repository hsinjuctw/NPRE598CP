import numpy as np
import Library.coord as coord

M  = -8e15 # dipole moment of Earth [Tm^3]
RE = 6.371e6

def dipoleEarthSph(r,theta,phi):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    #M  = -8.e15 # dipole moment of Earth [Tm^3]
    Br = 2.*M*np.cos(theta)/r**3
    Bt =    M*np.sin(theta)/r**3
    Bp = 0.
    return Br,Bt,Bp

def dipoleEarth(x,y,z):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    #M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(x**2+y**2+z**2)
    Bx = 3.*M*x*z/r**5
    By = 3.*M*y*z/r**5
    Bz = M*(3.*z*z-r*r)/r**5
    return Bx,By,Bz

# Y = 0 plane
def dipoleEarthX0Z(x,z):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    #M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(x**2+z**2)
    Bx = 3.*M*x*z/r**5
    #By = 3.*M*y*z/r**5
    Bz = M*(3.*z*z-r*r)/r**5
    return Bx,Bz

def dipoleEarth0YZ(y,z):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    #M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(y**2+z**2)
    #Bx = 3.*M*x*z/r**5
    By = 3.*M*y*z/r**5
    Bz = M*(3.*z*z-r*r)/r**5
    return By,Bz

def dipoleEarthXY0(x,y):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    #M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(x**2+y**2)
    Bx = 3.*M*x*z/r**5
    By = 3.*M*y*z/r**5
    #Bz = M*(3.*z*z-r*r)/r**5
    return Bx,By

def getEarthDipoleCSC(X,Y,Z):
    '''
    Get the dipole field (Bx, By, Bz, Bmag) of the Earth by transforming Cartesian to Spherical to Cartesian.
    '''
    Nx = len(X)
    Ny = len(Y)
    Nz = len(Z)
    Bmag = np.empty((Nx,Ny,Nz))
    Bx   = np.empty((Nx,Ny,Nz))
    By   = np.empty((Nx,Ny,Nz))
    Bz   = np.empty((Nx,Ny,Nz))
    for i in range(0,Nx):
        for j in range(0,Ny):
            for k in range(0,Nz):
                if abs(X[i])<RE and abs(Y[j])<RE and abs(Z[k])<RE:
                    Bmag[i,j,k] = np.nan
                    Bx[i,j,k]   = np.nan
                    By[i,j,k]   = np.nan
                    Bz[i,j,k]   = np.nan
                else:
                    r,theta,phi = coord.car2sph( X[i],Y[j],Z[k] )
                    Br,Bt,Bp    = dipoleEarthSph( r,theta,phi )
                    Bx[i,j,k],By[i,j,k],Bz[i,j,k] = coord.sph2carV(Br,Bt,Bp,r,theta,phi)
                    Bmag[i,j,k] = np.sqrt( Bx[i,j,k]**2 + By[i,j,k]**2 + Bz[i,j,k]**2 )
    return Bx,By,Bz,Bmag

def getEarthDipole(X,Y,Z):
    '''
    Get the dipole field (Bx, By, Bz, Bmag) of the Earth in Cartesian.
    '''
    Nx = len(X)
    Ny = len(Y)
    Nz = len(Z)
    Bmag = np.empty((Nx,Ny,Nz))
    Bx   = np.empty((Nx,Ny,Nz))
    By   = np.empty((Nx,Ny,Nz))
    Bz   = np.empty((Nx,Ny,Nz))    
    for i in range(0,Nx):
        for j in range(0,Ny):
            for k in range(0,Nz):
                if abs(X[i])<RE and abs(Y[j])<RE and abs(Z[k])<RE:
                    Bmag[i,j,k] = np.nan
                    Bx[i,j,k]   = np.nan
                    By[i,j,k]   = np.nan
                    Bz[i,j,k]   = np.nan
                else:
                    Bx[i,j,k],By[i,j,k],Bz[i,j,k] = dipoleEarth(X[i],Y[j],Z[k])
                    Bmag[i,j,k] = np.sqrt( Bx[i,j,k]**2 + By[i,j,k]**2 + Bz[i,j,k]**2 )
    return Bx,By,Bz,Bmag