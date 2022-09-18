import numpy as np

def sph2car(r,theta,phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

def car2sph(x,y,z):
    x2y2  = x**2 + y**2
    r     = np.sqrt(x2y2 + z**2)
    theta = np.arctan2(np.sqrt(x2y2),z)
    phi   = np.arctan2(y,x)
    return r,theta,phi

# https://stackoverflow.com/questions/4116658/faster-numpy-cartesian-to-spherical-coordinate-conversion
#def appendSpherical_np(xyz):
#    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
#    xy = xyz[:,0]**2 + xyz[:,1]**2
#    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
#    ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
#    #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
#    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
#    return ptsnew

def sph2carV(Ar,At,Ap,r,theta,phi):
    T = ([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi), np.cos(theta)],
         [np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)],
         [             -np.sin(phi),              np.cos(phi),             0])
    xyz = np.matmul(np.transpose(T),(Ar,At,Ap))
    Ax,Ay,Az = np.transpose(xyz)
    #Ax = np.sin(theta)*np.cos(phi)*Ar+np.cos(theta)*np.cos(phi)*At-np.sin(phi)*Ap
    #Ay = np.sin(theta)*np.sin(phi)*Ar+np.cos(theta)*np.sin(phi)*At+np.cos(phi)*Ap
    #Az =             np.cos(theta)*Ar            -np.sin(theta)*At
    #x, y, z  = sph2car(r,theta,phi)
    return Ax,Ay,Az#,x,y,z

def car2sphV(Ax,Ay,Az,x,y,z):
    r,theta,phi = car2sph(x,y,z)
    T = ([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi), np.cos(theta)],
         [np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)],
         [             -np.sin(phi),              np.cos(phi),             0])
    rtp = np.matmul(T,(Ax,Ay,Az))
    Ar,Atheta,Aphi     = np.transpose(rtp)
    return Ar,Atheta,Aphi#,r,theta,phi