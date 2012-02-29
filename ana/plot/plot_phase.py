import numpy as np
import matplotlib.cm as cm
import  matplotlib.pyplot as plt
file_index = 7
filename = '../dat/selc/selc_%3d.dat' % file_index
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
w2 = M_re**2 + M_im**2
phi = np.arctan2(M_im, M_re) * 180 / np.pi
w2[w2<1e-99] = 1e-99

x = M_re
y = M_im


figsize_x, figsize_y = 7, 6
params = {#'backend': 'ps',
          'axes.labelsize': 30,
          'text.fontsize': 20,
          'legend.fontsize': 20,
          'xtick.labelsize': 25,
          'ytick.labelsize': 25,
          'text.usetex': True,
          'figure.figsize': [figsize_x, figsize_y]}

plt.rcParams.update(params)

xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()

fig = plt.figure(figsize=(figsize_x, figsize_y))
fig.canvas.set_window_title(filename) 
#plt.subplots_adjust(hspace=0.5)
plt.subplots_adjust(bottom=0.15, right=1.0, top=0.9)
#plt.subplot(111)
ax = fig.add_subplot(111, aspect='equal')
plt.hexbin(x, y, gridsize=50, cmap=cm.jet, bins='log', C=w2)
plt.grid(True)
plt.xlim([-0.03, 0.03])
plt.ylim([-0.03, 0.03])
ax.set_xticks([-0.03, -0.015, 0.0, 0.015, 0.03])
ax.set_yticks([-0.03, -0.015, 0.0, 0.015, 0.03])
plt.xlabel(r'$\mathrm{Re}[M_p]$ (arb.)')
plt.ylabel(r'$\mathrm{Im}[M_p]$ (arb.)')
#plt.ylim([0.8, 4.2])

# #plt.axis([xmin, xmax, ymin, ymax])
# plt.subplot(212)
# plt.hexbin(x, y, gridsize=50, cmap=cm.jet, bins='log')
# #plt.axis([xmin, xmax, ymin, ymax])
# #plt.title("Hexagon binning with log color scale")
# #cb = plt.colorbar()
# #cb.set_label('log10(N)')
# plt.xlim([-0.02, 0.02])
# plt.ylim([-0.02, 0.02])
# #plt.ylim([0.8, 4.2])
plt.savefig( '../fig/corr_%3d.png' % file_index )

plt.show()
