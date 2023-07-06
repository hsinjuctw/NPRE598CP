import numpy as np
import matplotlib.pyplot as plt
import Library.bearth as bearth
import Library.particle as particle
import Library.coord as coord
from matplotlib.patches import Circle
from matplotlib import ticker

import plotly.express as px
import plotly.graph_objects as go

qe = 1.60217662e-19
me = 9.10938356e-31
mp = 1.6726219e-27
RE = 6.371e6 # Earth radius [m]

def main():
    ##################################
    ##### EARTH'S MAGNETIC FIELD #####
    ##################################
    # Nx = 101
    # Ny = 101
    # Nz = 101

    # X = np.linspace( -20.*RE, 20.*RE, Nx)
    # Y = np.linspace( -20.*RE, 20.*RE, Ny)
    # Z = np.linspace( -20.*RE, 20.*RE, Nz)
    
    # Bx,By,Bz,Bmag = bearth.getEarthDipole(X,Y,Z)

    ###############################
    ##### PARTICLE TRAJECTORY #####
    ###############################
    fcorrection = True # flag for frequency correction
    # charge [C], mass [kg]
    Q = qe
    M = mp*131.29 # Xe
    # M = mp*15.999 # O
    # M = mp*7
    # x0 = -10.*RE
    # y0 = 0.
    # z0 = 0.
    x0,y0,z0 = coord.sph2car(2.*RE,20./180.*np.pi,0.)
    # B field [nT]
    bx, by, bz = bearth.dipoleEarth(x0,y0,z0)
    b = np.sqrt(bx**2+by**2+bz**2)
    # cyclotron frequency [rad/s]
    # omega_c = np.abs(Q)*b/M
    # cyclotron period [s]
    # Tc = 2.*np.pi/omega_c

    Tev = 1.e4             # particle temperature in eV
    T = particle.eV2K(Tev) # particle temperature in K
    v = particle.E2v(Tev,M)
    # print(v)
    # initial velocity [m/s]
    vx0 = -bx*v/b
    vy0 = -by*v/b
    vz0 = -bz*v/b
    # vx0 = -1.e7
    # vy0 =  1.e7
    # vz0 =  1.e7
    # magnitude of v_perp
    # v0 = np.sqrt(vx0**2+vy0**2)
    # Larmor radius
    # r_L = v0/omega_c
    # initial position [m]
    #x0 += r_L
    # initial state vector
    BB0 = np.array((x0,y0,z0,vx0,vy0,vz0))
    # time grids [s]
    timetot = 3600
    # timetot = 200
    N_per_sec = 100
    time = np.linspace(0.,timetot,N_per_sec*timetot+1)
    dt = time[1]-time[0]
    # parameters
    params = np.array((dt, Q/M*dt/2.))
    # Boris integration
    BB = particle.dirBorisBunemann(time, BB0, params, fcorrection)
    
    ##############################
    ##### FIELD LINE TRACING #####
    ##############################
    # Lshell = 10.
    # MLT    = 12.
    dX  = .01*RE
    dY  = .01*RE
    dZ  = .01*RE

    FLx0,FLy0,FLz0 = bearth.dipoleFieldline3D( 10.*RE,0.,0.,dX,dY,dZ)
    FLx1,FLy1,FLz1 = bearth.dipoleFieldline3D(-10.*RE,0.,0.,dX,dY,dZ)
    FLx2,FLy2,FLz2 = bearth.dipoleFieldline3D( 5.*RE*np.sqrt(2), 5.*RE*np.sqrt(2),0.,dX,dY,dZ)
    FLx3,FLy3,FLz3 = bearth.dipoleFieldline3D(-5.*RE*np.sqrt(2), 5.*RE*np.sqrt(2),0.,dX,dY,dZ)
    FLx4,FLy4,FLz4 = bearth.dipoleFieldline3D( 5.*RE*np.sqrt(2),-5.*RE*np.sqrt(2),0.,dX,dY,dZ)
    FLx5,FLy5,FLz5 = bearth.dipoleFieldline3D(-5.*RE*np.sqrt(2),-5.*RE*np.sqrt(2),0.,dX,dY,dZ)
    FLx6,FLy6,FLz6 = bearth.dipoleFieldline3D(0., 10.*RE,0.,dX,dY,dZ)
    FLx7,FLy7,FLz7 = bearth.dipoleFieldline3D(0.,-10.*RE,0.,dX,dY,dZ)
    FLx = [FLx0,FLx1,FLx2,FLx3,FLx4,FLx5,FLx6,FLx7]
    FLy = [FLy0,FLy1,FLy2,FLy3,FLy4,FLy5,FLy6,FLy7]
    FLz = [FLz0,FLz1,FLz2,FLz3,FLz4,FLz5,FLz6,FLz7]

    #################
    ##### PLOTS #####
    #################
    fig = go.Figure(data=[go.Scatter3d(
        x=BB[:,0]/RE, y=BB[:,1]/RE, z=BB[:,2]/RE,
        mode='markers',
        marker=dict(
            size=1,
    #         color=tSec,             # set color to an array/list of desired values
            # colorscale='rdbu',   # choose a colorscale
            opacity=0.8,
            # showscale = False
        )
    )])
    for i in range(8):
        fig.add_trace(go.Scatter3d(
            x=FLx[i]/RE, y=FLy[i]/RE, z=FLz[i]/RE,
            mode='markers',
            marker=dict(
                size=1,
                color='black',
                opacity=0.8,
                # showscale = False
            )
            ))

    # RE = 1. # Earth radius
    phi = np.linspace(0, 2*np.pi)
    theta = np.linspace((-1*np.pi)/2, np.pi/2)
    phi, theta = np.meshgrid(phi, theta)
    # Earth's surface
    xe = np.cos(theta) * np.sin(phi)
    ye = np.cos(theta) * np.cos(phi)
    ze = np.sin(theta)

    fig.add_trace(go.Surface(
        x = xe, y = ye, z = ze,
        colorscale=[[0, 'dimgray'], [1, 'dimgray']],
        showscale = False))
    fig.update_layout(scene_aspectmode='data')
    fig.update_coloraxes(showscale=False)
    fig.show()
    # fig.write_image("Traj.png")

if __name__ == '__main__':
   main()