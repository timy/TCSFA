import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('dat/crf_map.dat')
data_re = data[:,0]
data_im = data[:,1]
data_hi = data[:,2]

#fig = plt.figure()
#ax = Axes3D(fig)

scatter(data_re, data_im, c=log10(data_hi));
plt.xlabel('Re')
plt.ylabel('Im')

grid(True)

plt.show()
