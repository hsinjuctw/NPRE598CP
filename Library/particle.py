import numpy as np
#import matplotlib.pyplot as plt
import Library.bearth as bearth
from scipy.special import erf

qe = 1.60217662e-19
me = 9.10938356e-31
mp = 1.6726219e-27
vc = 2.99792458e+8  #speed of light
kB = 1.380649e-23

def Efield(x,y,z):
    Ex = 0.
    Ey = 0.
    Ez = 0.
    return Ex,Ey,Ez

def frequency_correction(x):
    if x == 0.:
        alpha = 1.
    else:
        alpha = np.tan(x)/x
    return alpha

def lorentz(v):
    return 1/np.sqrt(1-v**2/vc**2)

def dirBorisBunemann( time, x0, params, fcorrection):
    '''
    This Boris-Bunemann routine calculates magnetic field in each iteration.
    '''
    dt     = params[0]
    qmdt2  = params[1] # (q/m)*(dt/2)
    N      = np.size(time)
    M      = np.size(x0)
    X      = np.zeros((N,M))
    # X      = np.zeros((N,M+4))
    X[0,:] = x0
    # X[0,0:6] = x0
    x      = X[0,0]
    y      = X[0,1]
    z      = X[0,2]
    vx     = X[0,3]
    vy     = X[0,4]
    vz     = X[0,5]
    # Bx,By,Bz = bearth.dipoleEarth(x,y,z)
    # X[0,6] = np.sqrt(Bx**2+By**2+Bz**2) # Bmag
    # vpar   = (vx*Bx+vy*By+vz*Bz)/X[0,6]
    # X[0,7] = vx-vpar*Bx/X[0,6]# vperpx
    # X[0,8] = vy-vpar*By/X[0,6]# vperpy
    # X[0,9] = vz-vpar*Bz/X[0,6]# vperpz
    for i in range(0,N): # 0 for velocity pushback (v[-1/2]), 1 to N for the N steps
        Ex,Ey,Ez = Efield(x,y,z)
        Bx,By,Bz = bearth.dipoleEarth(x,y,z)
        if i == 0:
            # frequency correction: replace Omega*dt/2 with tan(Omega*dt/2)
            if fcorrection == True:
                alpha = frequency_correction(qmdt2*np.sqrt(Bx*Bx+By*By+Bz*Bz))
                alpha_x = frequency_correction(qmdt2*Bx/2.)
                alpha_y = frequency_correction(qmdt2*By/2.)
                alpha_z = frequency_correction(qmdt2*Bz/2.)
            else:
                alpha = 1.
            # 1. half acceleration along E field: v = v[0]-(q/m)(dt/2/2)E (half timestep back)
            vx -= qmdt2 * Ex * alpha / 2.
            vy -= qmdt2 * Ey * alpha / 2.
            vz -= qmdt2 * Ez * alpha / 2.
            # 2. B field rotation (half timestep back)
            tx = -qmdt2 * Bx * alpha / 2.
            ty = -qmdt2 * By * alpha / 2.
            tz = -qmdt2 * Bz * alpha / 2.
        else:
            # frequency correction: replace Omega*dt/2 with tan(Omega*dt/2)
            if fcorrection == True:
                alpha = frequency_correction(qmdt2*np.sqrt(Bx*Bx+By*By+Bz*Bz))
                alpha_x = frequency_correction(qmdt2*Bx)
                alpha_y = frequency_correction(qmdt2*By)
                alpha_z = frequency_correction(qmdt2*Bz)
            else:
                alpha = 1.
            # 1. half acceleration along E field: v = v[n-1/2]+(q/m)(dt/2)E
            vx += qmdt2 * Ex * alpha
            vy += qmdt2 * Ey * alpha
            vz += qmdt2 * Ez * alpha
            # 2. B field rotation
            tx = qmdt2 * Bx * alpha
            ty = qmdt2 * By * alpha
            tz = qmdt2 * Bz * alpha
        tmagsq = tx**2 + ty**2 + tz**2
        sx = 2.*tx/(1+tmagsq)
        sy = 2.*ty/(1+tmagsq)
        sz = 2.*tz/(1+tmagsq)
        # |x  y  z |
        # |vx vy vz| = (vy*tz - vz*ty, vz*tx - vx*tz, vx*ty - vy*tx)
        # |tx ty tz|
        # vp = v + v x t
        vpx = vx + vy*tz - vz*ty
        vpy = vy + vz*tx - vx*tz
        vpz = vz + vx*ty - vy*tx
        # v = v + vp x s
        vx += vpy*sz - vpz*sy
        vy += vpz*sx - vpx*sz
        vz += vpx*sy - vpy*sx
        if i == 0:
            # 3. half acceleration along E field (half timestep back)
            vx -= qmdt2 * Ex * alpha / 2.
            vy -= qmdt2 * Ey * alpha / 2.
            vz -= qmdt2 * Ez * alpha / 2.
        else:
            # 3. half acceleration along E field
            vx += qmdt2 * Ex * alpha
            vy += qmdt2 * Ey * alpha
            vz += qmdt2 * Ez * alpha
            # 4. push position (no need to update for initial condition)
            x += vx*dt
            y += vy*dt
            z += vz*dt
        # store the coordinates back into X
        X[i,0] = x
        X[i,1] = y
        X[i,2] = z
        # v[i] in code is v[i-1/2] in theory
        X[i,3] = vx
        X[i,4] = vy
        X[i,5] = vz
        # X[i,6] = np.sqrt(Bx**2+By**2+Bz**2)
        # vpar   = (vx*Bx+vy*By+vz*Bz)/X[0,6]
        # X[i,7] = vx-vpar*Bx/X[i,6]# vperpx
        # X[i,8] = vy-vpar*By/X[i,6]# vperpy
        # X[i,9] = vz-vpar*Bz/X[i,6]# vperpz
    return X

def interp3D(X,Y,Z,Bx_grid,By_grid,Bz_grid,xp,yp,zp):
    '''
    a very versatile but super inefficient 3D interpolation routine for 3 3D vectors
    '''
    dX = X[1]-X[0]
    dY = Y[1]-Y[0]
    dZ = Z[1]-Z[0]
    # Obtain index along X, Y, Z coordinates
    i = int(np.absolute( np.floor( (xp-X[0])/dX )))
    j = int(np.absolute( np.floor( (yp-Y[0])/dY )))
    k = int(np.absolute( np.floor( (zp-Z[0])/dZ )))
    # Get the coordinates of the closest node along x,y
    xi = X[0] + i * dX
    yj = Y[0] + j * dY
    zk = Z[0] + k * dZ
    # Cell volumes
    # 0-7 are I-VIII octants
    A0 = (xi+dX-xp)*(yj+dY-yp)*(zk+dZ-zp)
    A1 = (xp-xi)   *(yj+dY-yp)*(zk+dZ-zp)
    A2 = (xp-xi)   *(yp-yj)   *(zk+dZ-zp)
    A3 = (xi+dX-xp)*(yp-yj)   *(zk+dZ-zp)
    A4 = (xi+dX-xp)*(yj+dY-yp)*(zp-zk)
    A5 = (xp-xi)   *(yj+dY-yp)*(zp-zk)
    A6 = (xp-xi)   *(yp-yj)   *(zp-zk)
    A7 = (xi+dX-xp)*(yp-yj)   *(zp-zk)
    At = dX * dY * dZ
    # Linear weights
    w0 = A0 / At
    w1 = A1 / At
    w2 = A2 / At
    w3 = A3 / At
    w4 = A4 / At
    w5 = A5 / At
    w6 = A6 / At
    w7 = A7 / At
    Bx = w0 * Bx_grid[ i ][ j ][ k ]
    By = w0 * By_grid[ i ][ j ][ k ]
    Bz = w0 * Bz_grid[ i ][ j ][ k ]
    if xi == X[-1]:
        if yj == Y[-1]:
            if zk != Z[-1]:
                Bx += w4 * Bx_grid[ i ][ j ][k+1]
                By += w4 * By_grid[ i ][ j ][k+1]
                Bz += w4 * Bz_grid[ i ][ j ][k+1]
        else:
            Bx += w3 * Bx_grid[ i ][j+1][ k ]
            By += w3 * By_grid[ i ][j+1][ k ]
            Bz += w3 * Bz_grid[ i ][j+1][ k ]
            if zk != Z[-1]:
                Bx += w4 * Bx_grid[ i ][ j ][k+1] + w7 * Bx_grid[ i ][j+1][k+1]
                By += w4 * By_grid[ i ][ j ][k+1] + w7 * By_grid[ i ][j+1][k+1]
                Bz += w4 * Bz_grid[ i ][ j ][k+1] + w7 * Bz_grid[ i ][j+1][k+1]
    else:
        Bx += w1 * Bx_grid[i+1][ j ][ k ]
        By += w1 * By_grid[i+1][ j ][ k ]
        Bz += w1 * Bz_grid[i+1][ j ][ k ]
        if yj == Y[-1]:
            if zk != Z[-1]:
                Bx += w4 * Bx_grid[ i ][ j ][k+1] + w5 * Bx_grid[i+1][ j ][k+1]
                By += w4 * By_grid[ i ][ j ][k+1] + w5 * By_grid[i+1][ j ][k+1]
                Bz += w4 * Bz_grid[ i ][ j ][k+1] + w5 * Bz_grid[i+1][ j ][k+1]
        else:
            Bx += w2 * Bx_grid[i+1][j+1][ k ] + w3 * Bx_grid[ i ][j+1][ k ]
            By += w2 * By_grid[i+1][j+1][ k ] + w3 * By_grid[ i ][j+1][ k ]
            Bz += w2 * Bz_grid[i+1][j+1][ k ] + w3 * Bz_grid[ i ][j+1][ k ]
            if zk != Z[-1]:
                Bx += w4 * Bx_grid[ i ][ j ][k+1] + w5 * Bx_grid[i+1][ j ][k+1] \
                    + w6 * Bx_grid[i+1][j+1][k+1] + w7 * Bx_grid[ i ][j+1][k+1]
                By += w4 * By_grid[ i ][ j ][k+1] + w5 * By_grid[i+1][ j ][k+1] \
                    + w6 * By_grid[i+1][j+1][k+1] + w7 * By_grid[ i ][j+1][k+1]
                Bz += w4 * Bz_grid[ i ][ j ][k+1] + w5 * Bz_grid[i+1][ j ][k+1] \
                    + w6 * Bz_grid[i+1][j+1][k+1] + w7 * Bz_grid[ i ][j+1][k+1]
        # full equation for xi != X[-1], yj != Y[-1], zk != Z[-1]
        # B = w0 * B_grid[  i,  j,  k] + w1 * B_grid[i+1,  j,  k] + w2 * B_grid[i+1,j+1,  k] + w3 * B_grid[  i,j+1,  k] \
        #   + w4 * B_grid[  i,  j,k+1] + w5 * B_grid[i+1,  j,k+1] + w6 * B_grid[i+1,j+1,k+1] + w7 * B_grid[  i,j+1,k+1]
    return Bx,By,Bz

def interp2D(X,Y,Bx_grid,By_grid,Bz_grid,xp,yp):
    '''
    a very versatile but super inefficient 2D interpolation routine for 3 2D vectors
    '''
    dX = X[1]-X[0]
    dY = Y[1]-Y[0]
    # Obtain index along X, Y coordinates
    i = int(np.absolute( np.floor( (xp-X[0])/dX )))
    j = int(np.absolute( np.floor( (yp-Y[0])/dY )))
    # Get the coordinates of the closest node along x,y
    xi = X[0] + i * dX
    yj = Y[0] + j * dY
    # Cell volumes
    # 0-7 are I-VIII octants
    A0 = (xi+dX-xp)*(yj+dY-yp)
    A1 = (xp-xi)   *(yj+dY-yp)
    A2 = (xp-xi)   *(yp-yj)
    A3 = (xi+dX-xp)*(yp-yj)
    At = dX * dY
    # Linear weights
    w0 = A0 / At
    w1 = A1 / At
    w2 = A2 / At
    w3 = A3 / At
    Bx = w0 * Bx_grid[i][j]
    By = w0 * By_grid[i][j]
    Bz = w0 * Bz_grid[i][j]
    if xi == X[-1]:
        if yj != Y[-1]:
            Bx += w3 * Bx_grid[i][j+1]
            By += w3 * By_grid[i][j+1]
            Bz += w3 * Bz_grid[i][j+1]
    else:
        Bx += w1 * Bx_grid[i+1][j]
        By += w1 * By_grid[i+1][j]
        Bz += w1 * Bz_grid[i+1][j]
        if yj != Y[-1]:
            Bx += w2 * Bx_grid[i+1][j+1] + w3 * Bx_grid[i][j+1]
            By += w2 * By_grid[i+1][j+1] + w3 * By_grid[i][j+1]
            Bz += w2 * Bz_grid[i+1][j+1] + w3 * Bz_grid[i][j+1]
    return Bx,By,Bz

def interp1D(X,Vx,xp):
    '''
    a very versatile but super inefficient 1D interpolation routine for 1D vector of non-uniform position spacing
    this uses twice as much computation resources than needed due to the symmetry of Maxwellian velocity distribution
    '''
    # Obtain index along X coordinates
    for i in range(0,len(X)):
        if xp <= X[i]:
            i = i-1
            break
    if i < 0:
        Vxi = X[0] # caps at extreme value
    elif X[i] == X[-1]:
        Vxi = Vx[-1] # caps at extreme value
    else:
        dX = X[i+1]-X[i]
        A0 = (X[i+1]-xp)
        A1 = (xp-X[i])
        # Linear weights
        w0 = A0 / dX
        w1 = A1 / dX
        Vxi = w0 * Vx[i] + w1 * Vx[i+1]
    return Vxi

def maxwellian1D(m,v,T):
    '''
    Returns probability density at 1D velocity v given mass m and temperature T.
    '''
    f = np.sqrt(m/(2.*np.pi*kB*T)) * np.exp(-m*v**2/(2.*kB*T))
    return f

# def findVelocity(percentile, m, T):
#     return -(2.*kB*T)*np.log(percentile/np.sqrt(m/(2.*np.pi*kB*T)))/m

# https://notebook.community/tommyogden/quantum-python-lectures/11_Monte-Carlo-Maxwell-Boltzmann-Distributions
def maxwellianCDF1D(m,v,T):
    return .5*erf(np.sqrt(m/2./kB/T)*v)+.5

def invMaxwellianCDF1D(m,T,Vx,Xp):
    '''
    Vx: a predefined vector of 1D velocities
    Xp: a vector of CDF values that needs interpolation
    '''
    cdf = maxwellianCDF1D(m,Vx,T)
    Vxi = np.empty(len(Xp)) # interpolated velocities
    for i in range(0,len(Xp)):
        Vxi[i] = interp1D(cdf,Vx,Xp[i])
    return Vxi

def eV2K(kBT):
    return kBT*qe/kB

def E2v(kBT,M):
    '''
    Returns speed (with relativistic correction) when given temperature in eV.
    '''
    E = kBT*qe
    v  = np.sqrt(-((1/(E/M/vc**2+1))**2-1)*vc**2)
    return v
