import numpy as np
import matplotlib.pyplot as plt
import Library.bearth as bearth
import Library.particle as particle

qe = 1.60217662e-19
me = 9.10938356e-31
#mp = 1.6726219e-27

def main():
    RE = 6.371e6 # Earth radius [m]
    ###############################
    ##### PARTICLE TRAJECTORY #####
    ###############################
    fcorrection = True # flag for frequency correction
    # charge [C], mass [kg]
    Q = -qe
    M = me
    x0 = -10.*RE
    y0 = 0.
    z0 = 0.
    # B field [nT]
    bx, by, bz = bearth.dipoleEarth(x0,y0,z0)
    b = np.sqrt(bx**2+by**2+bz**2)
    # cyclotron frequency [rad/s]
    omega_c = np.abs(Q)*b/M
    # cyclotron period [s]
    Tc = 2.*np.pi/omega_c
    # initial velocity [m/s]
    vx0 = -1.e7
    vy0 = 1.e7
    vz0 = 1.e7
    # magnitude of v_perp
    v0 = np.sqrt(vx0**2+vy0**2)
    # Larmor radius
    r_L = v0/omega_c
    # initial position [m]
    #x0 += r_L
    # initial state vector
    X0 = np.array((x0,y0,z0,vx0,vy0,vz0))
    # time grids [s]
    timetot = 1
    Ntime = 100
    time = np.linspace(0.,timetot,Ntime+1)
    dt = time[1]-time[0]
    # parameters
    params = np.array((dt, Q/M*dt/2.))
    # Boris integration
    X = particle.boris_bunemann(time, X0, params, fcorrection)
    
    # if fcorrection == True:
    #     plabel = 'Boris–Bunemann with frequency correction'
    # else:
    #     plabel = 'Boris–Bunemann without frequency correction'
    plt.figure(4)
    plt.plot( X[:,0]/RE, X[:,2]/RE, 'b-')#, label=plabel )
    plt.xlabel('$x$ [$R_E$]')
    plt.ylabel('$z$ [$R_E$]')
    #plt.axis('equal')
    plt.legend(loc=3)
    plt.savefig('boris_XZ.png')
    plt.figure(5)
    plt.plot( X[:,0]/RE, X[:,1]/RE, 'b-')#, label=plabel )
    plt.xlabel('$x$ [$R_E$]')
    plt.ylabel('$y$ [$R_E$]')
    plt.axis('equal')
    plt.legend(loc=3)
    plt.savefig('boris_XY.png')
    plt.figure(6)
    plt.plot( X[:,1]/RE, X[:,2]/RE, 'b-')#, label=plabel )
    plt.xlabel('$y$ [$R_E$]')
    plt.ylabel('$z$ [$R_E$]')
    plt.axis('equal')
    plt.legend(loc=3)
    plt.savefig('boris_YZ.png')
    plt.figure(7)
    plt.plot( time, X[:,0]/RE, 'r-', label='$x(t)$')
    plt.plot( time, X[:,1]/RE, 'g-', label='$y(t)$')
    plt.plot( time, X[:,2]/RE, 'b-', label='$z(t)$')
    plt.xlabel('$t$ [s]')
    plt.ylabel('location [$R_E$]')
    plt.legend(loc=3)
    plt.savefig('boris_time.png')
    plt.show()

if __name__ == '__main__':
   main()