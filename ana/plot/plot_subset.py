'''
This file is used to plot the correlation of subset

'''

from matplotlib import *
from pylab import *
import numpy as np

data = np.loadtxt('../../dat/dat_subset.dat')

px0, pz0 = data[:,0], data[:,1];
t0, ts_im = data[:,2], data[:,3];
x0, z0 = data[:,4], data[:,5];
px, pz = data[:,6], data[:,7];
L = data[:,8];
M_re, M_im = data[:,9], data[:,10];
n_pass_x, n_pass_z = data[:,11], data[:,12]


plot( t0, n_pass_x, 'o', markersize=4 )

plt.xlabel('t0')
plt.ylabel('L')
grid(True)

#xlim(-20, 30)
ylim(-0.5,2.5)
plt.show()
