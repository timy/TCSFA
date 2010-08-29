'''
This file is used to plot the pulse

'''

from matplotlib import *
from pylab import *
import numpy as np

data = np.loadtxt('dat/pulse.dat')
t = data[:,0];
Az = data[:,1];
Ax = data[:,2];
Ez = data[:,3];
Ex = data[:,4];

n_t = 200 # should agree with the n_t


plot( t, Ez, 'r.-' );
#plot( [0], x[0], 'bo')
plt.xlabel('z')
plt.ylabel('x')
grid(True)

plt.show()
