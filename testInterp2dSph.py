import numpy as np
from scipy.stats import maxwell
import matplotlib.pyplot as plt
import Library.particle as particle

# mean, var, skew, kurt = maxwell.stats(10,moments='mvsk')
# print(mean)
# vals = maxwell.ppf([0.001, 0.5, 0.999],20)
# print(vals)
# np.allclose([0.001, 0.5, 0.999], maxwell.cdf(vals))
# r = maxwell.rvs(20,size=1000)

# fig, ax = plt.subplots(1, 1)
# ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
# ax.legend(loc='best', frameon=False)
# plt.show()

me = 9.10938356e-31
mp = 1.6726219e-27
Tev = 0.025 # particle temperature in eV
T = particle.eV2K(Tev)
v = particle.E2v(Tev)
# percentile = .1
# print(particle.findVelocity(percentile,mp,T))
# print(particle.maxwellian1D(mp,0.,T))

# X = np.linspace(-10000.,10000.,2000)
# f = np.empty(len(X))
# for i in range(0,len(X)):
#     f[i] = particle.maxwellian1D(mp,X[i],T)

# fig = plt.figure(1)
# ax = fig.add_subplot(111)
# plt.plot( X,f, 'k-', label='Distribution' )
# plt.xlabel('$v$ [m/v]')
# plt.ylabel('[\%]')
# plt.legend(loc=3)
# plt.show()
Np = 10
Vx = np.linspace(-2.*v,2.*v,100)
Xp = np.random.rand(Np)
# print(particle.maxwellianCDF1D(mp,Vx,T))
print(particle.invMaxwellianCDF1D(mp,T,Vx,[.5]))
#print(particle.invMaxwellianCDF1D(mp,T,Vx,Xp))