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

    coor = 'S' # 'C' for Cartesian, 'S' for spherical

    if coor == 'C':
    # CARTESIAN CASE
        Nx = 21
        Ny = 21
        Nz = 21

        X = np.linspace( -20.*RE, 20.*RE, Nx)
        Y = np.linspace( -20.*RE, 20.*RE, Ny)
        Z = np.linspace( -20.*RE, 20.*RE, Nz)

        Bx,By,Bz,Bmag = bearth.getEarthDipole(X,Y,Z)

        refine  = 3.
        Xrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Nx-1)+1))
        Yrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Ny-1)+1))
        Zrefine = np.linspace( -20.*RE, 20.*RE, int(refine*(Nz-1)+1))

        # interpolated b field array initialization
        Bxi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
        Byi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))
        Bzi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))

        for i in range(0,int(refine*(Nx-1)+1)):
            for j in range(0,int(refine*(Ny-1)+1)):
                for k in range(0,int(refine*(Nz-1)+1)):
                    if abs(Xrefine[i])<RE and abs(Yrefine[j])<RE and abs(Zrefine[k])<RE:
                        Bxi[i,j,k] = np.nan
                        Byi[i,j,k] = np.nan
                        Bzi[i,j,k] = np.nan
                    else:
                        xp = Xrefine[i]
                        yp = Yrefine[j]
                        zp = Zrefine[k]
                        Bxi[i,j,k],Byi[i,j,k],Bzi[i,j,k] = particle.interp3D(X,Y,Z,Bx,By,Bz,xp,yp,zp)    
        plt.figure(1)
        XX,ZZ = np.meshgrid(X,Z)
        plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bx[:,10,:])
        plt.xlabel('$X$ [$R_E$]')
        plt.ylabel('$Z$ [$R_E$]')
        plt.axis('equal')
        plt.colorbar()
        plt.title('Bx in XZ plane [T] - Earth magnetic field')

        plt.figure(2)
        XXrefine,ZZrefine = np.meshgrid(Xrefine,Zrefine)
        plt.contourf(np.transpose(XXrefine)/RE,np.transpose(ZZrefine)/RE, Bxi[:,int(10*refine),:])
        plt.xlabel('$X$ [$R_E$]')
        plt.ylabel('$Z$ [$R_E$]')
        plt.axis('equal')
        plt.colorbar()
        plt.title('Interpolated Bx in XZ plane [T] - Earth magnetic field')
        plt.show()

    elif coor == 'S':
        # SPHERICAL CASE
        Nr = 21
        Nt = 21
        Np = 21

        R     = np.linspace( RE,   20.*RE, Nr )
        Theta = np.linspace( 0., 2.*np.pi, Nt )
        Phi   = np.linspace( 0., 2.*np.pi, Np )

        Br,Btheta,Bphi = bearth.getEarthDipoleSph(R,Theta)
        
        refine  = 3.
        Rrefine = np.linspace( RE,   20.*RE, int(refine*(Nr-1)+1))
        Trefine = np.linspace( 0., 2.*np.pi, int(refine*(Nt-1)+1))

        # interpolated b field array initialization
        Bri = np.zeros((int(refine*(Nr-1)+1), int(refine*(Nt-1)+1)))
        Bti = np.zeros((int(refine*(Nr-1)+1), int(refine*(Nt-1)+1)))
        Bpi = np.zeros((int(refine*(Nr-1)+1), int(refine*(Nt-1)+1)))

        for i in range(0,int(refine*(Nr-1)+1)):
            for j in range(0,int(refine*(Nt-1)+1)):
                if abs(Rrefine[i])<RE:
                    Bri[i,j] = np.nan
                    Bti[i,j] = np.nan
                    Bpi[i,j] = np.nan
                else:
                    rp = Rrefine[i]
                    tp = Trefine[j]
                    Bri[i,j],Bti[i,j],Bpi[i,j] = particle.interp2D(R,Theta,Br,Btheta,Bphi,rp,tp)

        # plt.figure(1)
        fig1, ax1 = plt.subplots(subplot_kw=dict(projection='polar'))
        RR,TT = np.meshgrid(R,Theta)
        # plt.contourf(np.transpose(RR)/RE,np.transpose(TT), Br)
        ax1.contourf(TT, RR/RE, np.log(abs(np.transpose(Br))))
        # plt.colorbar()
        plt.title('$B_r$ in $r\\theta$ plane [T] - Earth magnetic field')

        # plt.figure(2)
        fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'))
        RRrefine,TTrefine = np.meshgrid(Rrefine,Trefine)
        # plt.contourf(np.transpose(RRrefine)/RE,np.transpose(TTrefine), Bri)
        ax2.contourf(TTrefine, RRrefine/RE, np.log(abs(np.transpose(Bri))))
        # plt.colorbar()
        plt.title('Interpolated $B_r$ in $r\\theta$ plane [T] - Earth magnetic field')
        plt.show()

if __name__ == '__main__':
   main()