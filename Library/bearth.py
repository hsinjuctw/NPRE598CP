import numpy as np
import Library.coord as coord

M  = -8e15 # dipole moment of Earth [Tm^3]
RE = 6.371e6

def dipoleEarthSph(r,theta):
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
    Bz = -M/r**3
    return Bz

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
                    Br,Bt,Bp    = dipoleEarthSph( r,theta )
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

def getEarthDipoleSph(R,Theta):
    '''
    Get the dipole field (Br, Btheta, Bphi) of the Earth in Spherical.
    '''
    Nr = len(R)
    Nt = len(Theta)
    Bmag = np.empty((Nr,Nt))
    Br   = np.empty((Nr,Nt))
    Bt   = np.empty((Nr,Nt))
    Bp   = np.empty((Nr,Nt))
    for i in range(0,Nr):
        for j in range(0,Nt):
            if abs(R[i])<RE:
                Bmag[i,j] = np.nan
                Br[i,j]   = np.nan
                Bt[i,j]   = np.nan
                Bp[i,j]   = np.nan
            else:
                Br[i,j],Bt[i,j],Bp[i,j] = dipoleEarthSph(R[i],Theta[j])
    return Br,Bt,Bp

# NOT WORKING
def dipoleFieldline2DSph(r,theta,dr,rdtheta,nstep):
    R = [r]
    Theta = [theta]
    for i in range(0,nstep):
        Br,Bt,Bp = dipoleEarthSph(r,theta)
        theta += Bt*rdtheta/abs(Br)
        r     += Br*dr/abs(Br)
        R      = np.append(R,r)
        Theta  = np.append(Theta,theta)
        print(r,theta)
    return R,Theta

def dipoleFieldline3D(x,y,z,dx,dy,dz,nstep):
    X = [x]
    Y = [y]
    Z = [z]
    for i in range(0,nstep):
        Bx,By,Bz = dipoleEarth(x,y,z)
        Bmag     = np.sqrt(Bx**2+By**2+Bz**2)
        x       += dx*Bx/Bmag
        y       += dy*By/Bmag
        z       += dz*Bz/Bmag
        X        = np.append(X,x)
        Y        = np.append(Y,y)
        Z        = np.append(Z,z)
        print(x,y,z)
    return X,Y,Z

def dipoleFieldline2D(Lshell,MLT,dx,dz):
    x = Lshell*RE
    y = 0.
    z = 0.
    X = [x]
    Z = [z]
    # 1st quadrant
    while x**2+z**2 > RE**2:
        Bx,By,Bz = dipoleEarth(x,y,z)
        Bmag     = np.sqrt(Bx**2+By**2+Bz**2)
        x       += dx*Bx/Bmag
        # y       += dy*By/Bmag
        z       += dz*Bz/Bmag
        X        = np.append(X,x)
        # Y        = np.append(Y,y)
        Z        = np.append(Z,z)
    X = np.append(X[::-1],X)
    Z = np.append(Z[::-1],-Z)
    Y = np.sin((MLT-12)*15./180.*np.pi)*X
    X = np.cos((MLT-12)*15./180.*np.pi)*X
    return X,Y,Z




def blines(y,x, filament, current):
    X=y[0]
    Y=y[1]
    Z=y[2]
    direction=y[3]
    point = np.array([ [X], [Y], [Z] ])
    B     = biotsavart( filament, current, point )
    Bnorm = np.sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
    dY    = np.zeros(4)
    dY[0] = direction * B[0]/Bnorm
    dY[1] = direction * B[1]/Bnorm
    dY[2] = direction * B[2]/Bnorm
    dY[3] = 0.0
    return dY