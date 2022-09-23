import numpy as np

def sph2car(r,theta,phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

def car2sph(x,y,z):
    # https://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
    x2y2  = x**2 + y**2
    r     = np.sqrt(x2y2 + z**2)
    theta = np.arctan2(np.sqrt(x2y2),z)
    phi   = np.arctan2(y,x)
    return r,theta,phi

def sph2carV(Ar,At,Ap,r,theta,phi):
    T = ([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi), np.cos(theta)],
         [np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)],
         [             -np.sin(phi),              np.cos(phi),             0])
    xyz = np.matmul(np.transpose(T),(Ar,At,Ap))
    Ax,Ay,Az = np.transpose(xyz)
    return Ax,Ay,Az

def car2sphV(Ax,Ay,Az,x,y,z):
    r,theta,phi = car2sph(x,y,z)
    T = ([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi), np.cos(theta)],
         [np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)],
         [             -np.sin(phi),              np.cos(phi),             0])
    rtp = np.matmul(T,(Ax,Ay,Az))
    Ar,Atheta,Aphi     = np.transpose(rtp)
    return Ar,Atheta,Aphi