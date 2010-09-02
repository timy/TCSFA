'''
This file is used to plot the trajectory for a given p0

'''

from matplotlib import *
from pylab import *
import numpy as np


t0 = []
z0 = []
x0 = []
vz0 = []
vx0 = []
Mp_re = []
Mp_im = []
w = []
i_type = []

figure(figsize=(8,8), frameon=False)
ax = subplot(111, frameon=True)
data = np.loadtxt('../data/traj.dat')
t = data[:,0];
x = data[:,1];
vx = data[:,2];
z = data[:,3];
vz = data[:,4];
t0.append(t[0]);
z0.append(z[0]);
x0.append(x[0]);
vz0.append(vz[0]);
vx0.append(vx[0]);

plot( z, x );
scatter( z[0], x[0], c='r', s=40)



#plt.title( 'maximum of ionization probability: ' + str( amax(w)) )
#plt.xlabel('z')
#plt.ylabel('x')

plt.grid(True)
plt.savefig('/data/Document/Poster/EMMI_Workshop_2010/latex-poster/fig/traj.svg')
plt.show()

