import numpy as np

#M  = -8e15 # dipole moment of Earth [Tm^3]

def dipoleEarthSph(r,theta,phi):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    M  = -8.e15 # dipole moment of Earth [Tm^3]
    Br = 2.*M*np.cos(theta)/r**3
    Bt =    M*np.sin(theta)/r**3
    Bp = 0.
    return Br,Bt,Bp

def dipoleEarth(x,y,z):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(x**2+y**2+z**2)
    Bx = 3.*M*x*z/r**5
    By = 3.*M*y*z/r**5
    Bz = M*(3.*z*z-r*r)/r**5
    return Bx,By,Bz

# Y = 0 plane
def dipoleEarthX0Z(x,z):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(x**2+z**2)
    Bx = 3.*M*x*z/r**5
    #By = 3.*M*y*z/r**5
    Bz = M*(3.*z*z-r*r)/r**5
    return Bx,Bz

def dipoleEarth0YZ(y,z):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(y**2+z**2)
    #Bx = 3.*M*x*z/r**5
    By = 3.*M*y*z/r**5
    Bz = M*(3.*z*z-r*r)/r**5
    return By,Bz

def dipoleEarthXY0(x,y):
    # reference: eq.1-3 https://ccmc.gsfc.nasa.gov/static/files/Dipole.pdf
    M  = -8.e15 # dipole moment of Earth [Tm^3]
    r = np.sqrt(x**2+y**2)
    Bx = 3.*M*x*z/r**5
    By = 3.*M*y*z/r**5
    #Bz = M*(3.*z*z-r*r)/r**5
    return Bx,By