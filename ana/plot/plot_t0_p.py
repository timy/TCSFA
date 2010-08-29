'''
This file is used to plot the t0 vs. p_//

'''

from matplotlib import *
from pylab import *
from numpy import *

data = loadtxt("../dat/series 1.dat");
#data = loadtxt("../dat/select 1.dat");
n_traj = len(data)
data_px_0 = data[:,0]
data_pz_0 = data[:,1]
data_ts_re = data[:,2]
data_ts_im = data[:,3]
data_x_0 = data[:,4]
data_z_0 = data[:,5]
data_px = data[:,6]
data_pz = data[:,7]
data_M_re = data[:,8]
data_M_im = data[:,9]
data_n_pass_x = data[:,10]
data_n_pass_z = data[:,11]
data_err = data[:,12]
data_type = data[:,13]

print 'number of trajectories:', n_traj

w = []
for i_traj in range(n_traj):
    w.append(data_M_re[i_traj]**2 + data_M_im[i_traj]**2)

ratio = 400 / ( amax(log10(w)) - amin(log10(w)) )
w_size = ceil( ratio * ( log10( w ) - amin( log10(w) ) + 1 ) )

scatter( data_pz, data_px_0, c=data_type )
#scatter( sqrt(data_pz**2+data_px**2), data_ts_re, c=data_type )#, s=w_size, c=data_type )
#scatter( data_pz_0, log10(data_M_re**2+data_M_im**2), c=data_type ) # s=w_size)
#xlim([410,420])
#ylim([10.2,10.8])
#plt.title( 't0 v.s. p_parallel ' + str( amax(w)) )
plt.xlabel('pz')
plt.ylabel('p0_x')
grid(True)

plt.show()
