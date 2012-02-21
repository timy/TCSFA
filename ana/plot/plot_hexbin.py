import numpy as np
import matplotlib.cm as cm
import  matplotlib.pyplot as plt

filename = '../dat/selc_  4.dat'
data = np.loadtxt( filename )
px_0 = data[:, 0]
pz_0 = data[:, 1]
ts_re = data[:, 2]
ts_im = data[:, 3]
x0 = data[:, 4]
z0 = data[:, 5]
px = data[:, 6]
pz = data[:, 7]
L = data[:, 8]
M_re = data[:, 9]
M_im = data[:, 10]
iType = data[:, 12]
w = M_re**2 + M_im**2
w[w<1e-99] = 1e-99

x = z0
y = iType

xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()

plt.subplots_adjust(hspace=0.5)
plt.subplot(111)
plt.hexbin(x, y, C=w, gridsize=50, bins='log', cmap=cm.jet)
plt.axis([xmin, xmax, ymin, ymax])
plt.title("Hexagon binning with log color scale")
cb = plt.colorbar()
cb.set_label('log10(N)')

plt.show()
