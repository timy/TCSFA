import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('../data/im_integrand.dat')
ts_im = data[:,0]
int_re = data[:,1]
int_im = data[:,2]

fig = plt.figure()

plt.plot(ts_im, int_im, 'o-.')
plt.xlabel('ts_im')
plt.ylabel('int_im')

grid(True)

plt.show()
