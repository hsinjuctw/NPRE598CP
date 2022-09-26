import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import Library.bearth as bearth
import Library.coord as coord
import Library.particle as particle
from matplotlib import ticker
import time

def main():
    # record start time
    start = time.time()
    ##################################
    ##### EARTH'S MAGNETIC FIELD #####
    ##################################
    RE = 6.371e6 # Earth radius [m]

    coor = 'C' # 'C' for Cartesian, 'S' for spherical

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
        Bmi     = np.zeros((int(refine*(Nx-1)+1), int(refine*(Ny-1)+1), int(refine*(Nz-1)+1)))

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
                        Bmi[i,j,k] = np.sqrt(Bxi[i,j,k]**2+Byi[i,j,k]**2+Bzi[i,j,k]**2)
        
        # record end time
        end = time.time()
        print('The time of execution of above program is :',(end-start) * 10**3, ' ms')

        cmap = 'PRGn'
        
        plt.rc('legend',fontsize=12)
        plt.rc('axes',labelsize=12)
        plt.rc('xtick',labelsize=12)
        plt.rc('ytick',labelsize=12)

        fig, ax1 = plt.subplots(figsize=(4.8, 4.8))
        XX,ZZ = np.meshgrid(X,Z)
        plt.contourf(np.transpose(XX)/RE,np.transpose(ZZ)/RE, Bmag[:,10,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
        plt.xlabel('$X$ [$R_E$]')
        plt.ylabel('$Z$ [$R_E$]')
        # cbar = plt.colorbar(format='%.0e')
        # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
        ax1.set_aspect('equal', 'box')
        plt.tight_layout()
        # plt.savefig('xz_fine.png',dpi=150)
        plt.savefig('xz_coarse.png',dpi=150)
        # plt.title('B in XZ plane [T] - Earth magnetic field')

        fig, ax2 = plt.subplots(figsize=(4.8, 4.8))
        XXrefine,ZZrefine = np.meshgrid(Xrefine,Zrefine)
        plt.contourf(np.transpose(XXrefine)/RE,np.transpose(ZZrefine)/RE, Bmi[:,30,:],locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
        plt.xlabel('$X$ [$R_E$]')
        plt.ylabel('$Z$ [$R_E$]')
        # cbar = plt.colorbar(format='%.0e')
        # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
        ax2.set_aspect('equal', 'box')
        plt.tight_layout()
        plt.savefig('xz_refine.png',dpi=150)
        # plt.title('Interpolated B in XZ plane [T] - Earth magnetic field')
        plt.show()

    elif coor == 'S':
        # SPHERICAL CASE
        Nr = 61
        Nt = 61
        Np = 61

        R     = np.linspace( RE,   20.*RE, Nr )
        Theta = np.linspace( 0., 2.*np.pi, Nt )
        Phi   = np.linspace( 0., 2.*np.pi, Np )

        Br,Btheta,Bphi = bearth.getEarthDipoleSph(R,Theta)

        # refine  = 3.
        # Rrefine = np.linspace( RE,   20.*RE, int(refine*(Nr-1)+1))
        # Trefine = np.linspace( 0., 2.*np.pi, int(refine*(Nt-1)+1))

        # # interpolated b field array initialization
        # Bri = np.zeros((int(refine*(Nr-1)+1), int(refine*(Nt-1)+1)))
        # Bti = np.zeros((int(refine*(Nr-1)+1), int(refine*(Nt-1)+1)))
        # Bpi = np.zeros((int(refine*(Nr-1)+1), int(refine*(Nt-1)+1)))

        # for i in range(0,int(refine*(Nr-1)+1)):
        #     for j in range(0,int(refine*(Nt-1)+1)):
        #         if abs(Rrefine[i])<RE:
        #             Bri[i,j] = np.nan
        #             Bti[i,j] = np.nan
        #             Bpi[i,j] = np.nan
        #         else:
        #             rp = Rrefine[i]
        #             tp = Trefine[j]
        #             Bri[i,j],Bti[i,j],Bpi[i,j] = particle.interp2D(R,Theta,Br,Btheta,Bphi,rp,tp)

        # record end time
        end = time.time()
        print('The time of execution of above program is :',(end-start) * 10**3, ' ms')

        # plt.figure(1)
        cmap = 'PRGn'
        
        plt.rc('legend',fontsize=12)
        plt.rc('axes',labelsize=12)
        plt.rc('xtick',labelsize=12)
        plt.rc('ytick',labelsize=12)
        fig1, ax1 = plt.subplots(subplot_kw=dict(projection='polar'),figsize=(4.8, 4.8))
        RR,TT = np.meshgrid(R,Theta)
        # plt.plot( FLt,FLr/RE, 'k-' )
        # plt.contourf(np.transpose(RR)/RE,np.transpose(TT), Br)
        cf1 = ax1.contourf(TT, RR/RE, np.transpose(np.sqrt(Br**2+Btheta**2+Bphi**2)),locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
        ax1.set_theta_zero_location('N')
        cbar = fig1.colorbar(cf1, format='%.0e')
        cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
        # plt.title('$B$ in $r\\theta$ plane [T] - Earth magnetic field')
        plt.tight_layout()
        plt.savefig('rtheta_fine.png',dpi=150)

        # # plt.figure(2)
        # fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'),figsize=(4.8, 4.8))
        # RRrefine,TTrefine = np.meshgrid(Rrefine,Trefine)
        # # plt.contourf(np.transpose(RRrefine)/RE,np.transpose(TTrefine), Bri)
        # cf2 = ax2.contourf(TTrefine, RRrefine/RE, np.transpose(np.sqrt(Bri**2+Bti**2+Bpi**2)),locator=ticker.LogLocator(),alpha=.5,cmap=cmap)
        # ax2.set_theta_zero_location('N')
        # cbar = fig2.colorbar(cf2, format='%.0e')
        # cbar.set_label('Magnetic field strength (T)', rotation=270, labelpad=15)
        # # plt.title('Interpolated $|B|$ in $r\\theta$ plane [T] - Earth magnetic field')
        # plt.tight_layout()
        # plt.savefig('rtheta_refined.png',dpi=150)
        plt.show()

if __name__ == '__main__':
   main()