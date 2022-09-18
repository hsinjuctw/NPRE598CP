from os import wait4
import numpy as np
#import matplotlib.pyplot as plt
import Library.bearth as bearth

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

def boris_bunemann( time, x0, params, fcorrection):#, Bfield, Efield):# assume static B and E fields
    dt     = params[0]
    qmdt2  = params[1] # (q/m)*(dt/2)
    N      = np.size(time)
    M      = np.size(x0)
    X      = np.zeros((N,M))
    X[0,:] = x0
    x      = X[0,0]
    y      = X[0,1]
    z      = X[0,2]
    vx     = X[0,3]
    vy     = X[0,4]
    vz     = X[0,5]
    for i in range(0,N): # 0 for velocity pushback (v[-1/2]), 1 to N for the N steps
        Ex,Ey,Ez = Efield(x,y,z)
        Bx,By,Bz = bearth.dipoleEarth(x,y,z)
        # frequency correction: replace Omega*dt/2 with tan(Omega*dt/2)
        if fcorrection == True:
            #alpha = frequency_correction(qmdt2*np.sqrt(Bx*Bx+By*By+Bz*Bz))
            alpha_x = frequency_correction(qmdt2*Bx)
            alpha_y = frequency_correction(qmdt2*By)
            alpha_z = frequency_correction(qmdt2*Bz)
        else:
            #alpha = 1.
            alpha_x = 1.
            alpha_y = 1.
            alpha_z = 1.
        if i == 0:
            # 1. half acceleration along E field: v = v[0]-(q/m)(dt/2/2)E (half timestep back)
            vx -= qmdt2 * Ex * alpha_x / 2.
            vy -= qmdt2 * Ey * alpha_y / 2.
            vz -= qmdt2 * Ez * alpha_z / 2.
            # 2. B field rotation (half timestep back)
            tx = -qmdt2 * Bx * alpha_x / 2.
            ty = -qmdt2 * By * alpha_y / 2.
            tz = -qmdt2 * Bz * alpha_z / 2.
        else:
            # 1. half acceleration along E field: v = v[n-1/2]+(q/m)(dt/2)E
            vx += qmdt2 * Ex * alpha_x
            vy += qmdt2 * Ey * alpha_y
            vz += qmdt2 * Ez * alpha_z
            # 2. B field rotation
            tx = qmdt2 * Bx * alpha_x
            ty = qmdt2 * By * alpha_y
            tz = qmdt2 * Bz * alpha_z
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
            vx -= qmdt2 * Ex * alpha_x / 2.
            vy -= qmdt2 * Ey * alpha_y / 2.
            vz -= qmdt2 * Ez * alpha_z / 2.
        else:
            # 3. half acceleration along E field
            vx += qmdt2 * Ex * alpha_x
            vy += qmdt2 * Ey * alpha_y
            vz += qmdt2 * Ez * alpha_z
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
    return X

def boris_bunemann_interp( time, x0, params, fcorrection, Bfield):
    dt     = params[0]
    qmdt2  = params[1] # (q/m)*(dt/2)
    N      = np.size(time)
    M      = np.size(x0)
    X      = np.zeros((N,M))
    X[0,:] = x0
    x      = X[0,0]
    y      = X[0,1]
    z      = X[0,2]
    vx     = X[0,3]
    vy     = X[0,4]
    vz     = X[0,5]
    for i in range(0,N): # 0 for velocity pushback (v[-1/2]), 1 to N for the N steps
        Ex,Ey,Ez = Efield(x,y,z)
        Bx,By,Bz = bearth.dipoleEarth(x,y,z)
        # frequency correction: replace Omega*dt/2 with tan(Omega*dt/2)
        if fcorrection == True:
            #alpha = frequency_correction(qmdt2*np.sqrt(Bx*Bx+By*By+Bz*Bz))
            alpha_x = frequency_correction(qmdt2*Bx)
            alpha_y = frequency_correction(qmdt2*By)
            alpha_z = frequency_correction(qmdt2*Bz)
        else:
            #alpha = 1.
            alpha_x = 1.
            alpha_y = 1.
            alpha_z = 1.
        if i == 0:
            # 1. half acceleration along E field: v = v[0]-(q/m)(dt/2/2)E (half timestep back)
            vx -= qmdt2 * Ex * alpha_x / 2.
            vy -= qmdt2 * Ey * alpha_y / 2.
            vz -= qmdt2 * Ez * alpha_z / 2.
            # 2. B field rotation (half timestep back)
            tx = -qmdt2 * Bx * alpha_x / 2.
            ty = -qmdt2 * By * alpha_y / 2.
            tz = -qmdt2 * Bz * alpha_z / 2.
        else:
            # 1. half acceleration along E field: v = v[n-1/2]+(q/m)(dt/2)E
            vx += qmdt2 * Ex * alpha_x
            vy += qmdt2 * Ey * alpha_y
            vz += qmdt2 * Ez * alpha_z
            # 2. B field rotation
            tx = qmdt2 * Bx * alpha_x
            ty = qmdt2 * By * alpha_y
            tz = qmdt2 * Bz * alpha_z
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
            vx -= qmdt2 * Ex * alpha_x / 2.
            vy -= qmdt2 * Ey * alpha_y / 2.
            vz -= qmdt2 * Ez * alpha_z / 2.
        else:
            # 3. half acceleration along E field
            vx += qmdt2 * Ex * alpha_x
            vy += qmdt2 * Ey * alpha_y
            vz += qmdt2 * Ez * alpha_z
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
    return X

# this currently doesn't work
def interp3Dalt(X,Y,Z,Bx_grid,By_grid,Bz_grid,xp,yp,zp):
    '''
    a very versatile but super inefficient 3D interpolation routine for 3D vectors
    '''
    RE = 6.371e6
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
    print('xi,yj,zk',xi/RE,yj/RE,zk/RE)
    # interpolate along X
    # 0                  1
    # xi----------xp-----xi+dX
    #     xp-xi    xi+dX-xp
    w0 = (xi + dX - xp)/dX
    w1 = (xp - xi)/dX
    #print('w0,w1',w0,w1)
    if xi == X[-1]:
        # on line y = yj   , z = zk
        Bx00 = Bx_grid[i][ j ][ k ]
        By00 = By_grid[i][ j ][ k ]
        Bz00 = Bz_grid[i][ j ][ k ]
        # on line y = yj+dY, z = zk
        Bx10 = Bx_grid[i][j+1][ k ]
        By10 = By_grid[i][j+1][ k ]
        Bz10 = Bz_grid[i][j+1][ k ]
        # on line y = yj   , z = zk+dZ
        Bx01 = Bx_grid[i][ j ][k+1]
        By01 = By_grid[i][ j ][k+1]
        Bz01 = Bz_grid[i][ j ][k+1]
        # on line y = yj+dY, z = zk+dZ
        Bx11 = Bx_grid[i][j+1][k+1]
        By11 = By_grid[i][j+1][k+1]
        Bz11 = Bz_grid[i][j+1][k+1]
    else:
        # on line y = yj   , z = zk
        Bx00 = w0*Bx_grid[i][ j ][ k ]+w1*Bx_grid[i+1][ j ][ k ]
        By00 = w0*By_grid[i][ j ][ k ]+w1*By_grid[i+1][ j ][ k ]
        Bz00 = w0*Bz_grid[i][ j ][ k ]+w1*Bz_grid[i+1][ j ][ k ]
        # on line y = yj+dY, z = zk
        Bx10 = w0*Bx_grid[i][j+1][ k ]+w1*Bx_grid[i+1][j+1][ k ]
        By10 = w0*By_grid[i][j+1][ k ]+w1*By_grid[i+1][j+1][ k ]
        Bz10 = w0*Bz_grid[i][j+1][ k ]+w1*Bz_grid[i+1][j+1][ k ]
        # on line y = yj   , z = zk+dZ
        Bx01 = w0*Bx_grid[i][ j ][k+1]+w1*Bx_grid[i+1][ j ][k+1]
        By01 = w0*By_grid[i][ j ][k+1]+w1*By_grid[i+1][ j ][k+1]
        Bz01 = w0*Bz_grid[i][ j ][k+1]+w1*Bz_grid[i+1][ j ][k+1]
        # on line y = yj+dY, z = zk+dZ
        Bx11 = w0*Bx_grid[i][j+1][k+1]+w1*Bx_grid[i+1][j+1][k+1]
        By11 = w0*By_grid[i][j+1][k+1]+w1*By_grid[i+1][j+1][k+1]
        Bz11 = w0*Bz_grid[i][j+1][k+1]+w1*Bz_grid[i+1][j+1][k+1]
    # interpolate along Y
    # 0                  1
    # yj----------yp-----yj+dY
    #     yp-yj    yj+dY-yp
    w2 = (yj + dY - yp)/dY
    w3 = (yp - yj)/dY
    #print('w2,w3',w2,w3)
    # on line z = zk
    Bx0 = w2*Bx00 + w3*Bx10
    By0 = w2*By00 + w3*By10
    Bz0 = w2*Bz00 + w3*Bz10
    # on line z = zk+dZ
    Bx1 = w2*Bx01 + w3*Bx11
    By1 = w2*By01 + w3*By11
    Bz1 = w2*Bz01 + w3*Bz11
    # interpolate along Z
    # 0                  1
    # zk----------zp-----zk+dZ
    #     zp-zk    zk+dZ-zp
    w4 = (zk + dZ - zp)/dZ
    w5 = (zp - zk)/dZ
    #print('w4,w5',w4,w5)
    Bx = w4*Bx0 + w5*Bx1
    By = w4*By0 + w5*By1
    Bz = w4*Bz0 + w5*Bz1
    return Bx,By,Bz

def interp3D(X,Y,Z,Bx_grid,By_grid,Bz_grid,xp,yp,zp):
    '''
    a very versatile but super inefficient 3D interpolation routine for 3D vectors
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
    #print(w0+w1+w2+w3+w4+w5+w6+w7)
    if xi == X[-1]:
        if yj == Y[-1]:
            Bx = w0 * Bx_grid[ i ][ j ][ k ] + w4 * Bx_grid[ i ][ j ][k+1]
            By = w0 * By_grid[ i ][ j ][ k ] + w4 * By_grid[ i ][ j ][k+1]
            Bz = w0 * Bz_grid[ i ][ j ][ k ] + w4 * Bz_grid[ i ][ j ][k+1]
        elif zk == Z[-1]:
            Bx = w0 * Bx_grid[ i ][ j ][ k ] + w3 * Bx_grid[ i ][j+1][ k ]
            By = w0 * By_grid[ i ][ j ][ k ] + w3 * By_grid[ i ][j+1][ k ]
            Bz = w0 * Bz_grid[ i ][ j ][ k ] + w3 * Bz_grid[ i ][j+1][ k ]
        else:
            Bx = w0 * Bx_grid[ i ][ j ][ k ] + w3 * Bx_grid[ i ][j+1][ k ] \
               + w4 * Bx_grid[ i ][ j ][k+1] + w7 * Bx_grid[ i ][j+1][k+1]
            By = w0 * By_grid[ i ][ j ][ k ] + w3 * By_grid[ i ][j+1][ k ] \
               + w4 * By_grid[ i ][ j ][k+1] + w7 * By_grid[ i ][j+1][k+1]
            Bz = w0 * Bz_grid[ i ][ j ][ k ] + w3 * Bz_grid[ i ][j+1][ k ] \
               + w4 * Bz_grid[ i ][ j ][k+1] + w7 * Bz_grid[ i ][j+1][k+1]
    elif yj == Y[-1]:
        if zk == Z[-1]:
            Bx = w0 * Bx_grid[ i ][ j ][ k ] + w1 * Bx_grid[i+1][ j ][ k ]
            By = w0 * By_grid[ i ][ j ][ k ] + w1 * By_grid[i+1][ j ][ k ]
            Bz = w0 * Bz_grid[ i ][ j ][ k ] + w1 * Bz_grid[i+1][ j ][ k ]
        else:
            Bx = w0 * Bx_grid[ i ][ j ][ k ] + w1 * Bx_grid[i+1][ j ][ k ] \
               + w4 * Bx_grid[ i ][ j ][k+1] + w5 * Bx_grid[i+1][ j ][k+1]
            By = w0 * By_grid[ i ][ j ][ k ] + w1 * By_grid[i+1][ j ][ k ] \
               + w4 * By_grid[ i ][ j ][k+1] + w5 * By_grid[i+1][ j ][k+1]
            Bz = w0 * Bz_grid[ i ][ j ][ k ] + w1 * Bz_grid[i+1][ j ][ k ] \
               + w4 * Bz_grid[ i ][ j ][k+1] + w5 * Bz_grid[i+1][ j ][k+1]
    elif zk == Z[-1]:
        Bx = w0 * Bx_grid[ i ][ j ][ k ] + w1 * Bx_grid[i+1][ j ][ k ] \
           + w2 * Bx_grid[i+1][j+1][ k ] + w3 * Bx_grid[ i ][j+1][ k ]
        By = w0 * By_grid[ i ][ j ][ k ] + w1 * By_grid[i+1][ j ][ k ] \
           + w2 * By_grid[i+1][j+1][ k ] + w3 * By_grid[ i ][j+1][ k ]
        Bz = w0 * Bz_grid[ i ][ j ][ k ] + w1 * Bz_grid[i+1][ j ][ k ] \
           + w2 * Bz_grid[i+1][j+1][ k ] + w3 * Bz_grid[ i ][j+1][ k ]
    else:
        Bx = w0 * Bx_grid[ i ][ j ][ k ] + w1 * Bx_grid[i+1][ j ][ k ] \
           + w2 * Bx_grid[i+1][j+1][ k ] + w3 * Bx_grid[ i ][j+1][ k ] \
           + w4 * Bx_grid[ i ][ j ][k+1] + w5 * Bx_grid[i+1][ j ][k+1] \
           + w6 * Bx_grid[i+1][j+1][k+1] + w7 * Bx_grid[ i ][j+1][k+1]
        By = w0 * By_grid[ i ][ j ][ k ] + w1 * By_grid[i+1][ j ][ k ] \
           + w2 * By_grid[i+1][j+1][ k ] + w3 * By_grid[ i ][j+1][ k ] \
           + w4 * By_grid[ i ][ j ][k+1] + w5 * By_grid[i+1][ j ][k+1] \
           + w6 * By_grid[i+1][j+1][k+1] + w7 * By_grid[ i ][j+1][k+1]
        Bz = w0 * Bz_grid[ i ][ j ][ k ] + w1 * Bz_grid[i+1][ j ][ k ] \
           + w2 * Bz_grid[i+1][j+1][ k ] + w3 * Bz_grid[ i ][j+1][ k ] \
           + w4 * Bz_grid[ i ][ j ][k+1] + w5 * Bz_grid[i+1][ j ][k+1] \
           + w6 * Bz_grid[i+1][j+1][k+1] + w7 * Bz_grid[ i ][j+1][k+1]
    return Bx,By,Bz

# this currently indexes out of bound
def interp2D(X,Y,Bx_grid,By_grid,Bz_grid,xp,yp):
    '''
    a very versatile but super inefficient 2D interpolation routine for 3D vectors
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
    Bx = w0 * Bx_grid[i][j] + w1 * Bx_grid[i+1][j] + w2 * Bx_grid[i+1][j+1] + w3 * Bx_grid[i][j+1]
    By = w0 * By_grid[i][j] + w1 * By_grid[i+1][j] + w2 * By_grid[i+1][j+1] + w3 * By_grid[i][j+1]
    Bz = w0 * Bz_grid[i][j] + w1 * Bz_grid[i+1][j] + w2 * Bz_grid[i+1][j+1] + w3 * Bz_grid[i][j+1]
    return Bx,By,Bz


def maxwellian1D(m,v,T):
    f = np.sqrt(m/(2.*np.pi*kB*T)) * np.exp(-m*v**2/(2.*kB*T))
    return f
'''
def main():
    kBT = 0.025
    T  = kBT*qe/1.38e-23
    print(T, ' K')
    print(kBT*qe,' J')
    E = kBT*qe
    v  = np.sqrt(-((1/(E/mp/vc**2+1))**2-1)*vc**2)
    print(v,' m/s')
    print(maxwellian1D(mp,v,T))

if __name__ == '__main__':
   main()
'''