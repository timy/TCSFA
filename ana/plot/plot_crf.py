import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *

data = np.loadtxt('../data/crf_map.dat')
data_re = data[:,0]
data_im = data[:,1]
data_hi = data[:,2]

n_grid_x = 160 # should agree with the n_grid in file 'find_root.f90'
n_grid_y = 80


#transfer to 2D matrix
arry_re = data_re.reshape(n_grid_y, n_grid_x)
arry_im = data_im.reshape(n_grid_y, n_grid_x)
arry_hi = log10( data_hi.reshape(n_grid_y, n_grid_x) )

print('re_min', np.min(arry_re))
print('re_max', np.max(arry_re))
print('im_min', np.min(arry_im))
print('im_max', np.max(arry_im))
print('w_min', np.min(arry_hi))
print('w_max',np.max(arry_hi))

fig = plt.figure(figsize=(16,6), frameon=False)
ax = fig.add_subplot(111)
levels = arange(-3, 2, 0.2)
cset1 = plt.contourf(arry_re, arry_im, arry_hi, levels,
                        cmap=cm.get_cmap('jet', len(levels)-1)
                        )

#cset1 = plt.contourf(arry_re, arry_im, arry_hi, 100, cmap=cm.jet, locator=ticker.LogLocator())



cs = plt.contour(arry_re, arry_im, arry_hi, levels,
                 linewidths=np.arange(0, 3, .5),
                 colors=('r', 'green', 'blue', (1,1,0), '#afeeee', '0.5')
                 )

colorbar(cset1)
#plt.clabel(cs, fontsize=9, inline=1)
plt.title('Searching the root on the complex plane')




data = np.loadtxt('../data/crf_track.dat')
re_pt = data[:,0]
im_pt = data[:,1]
dr_pt = data[:,2]


scatter(re_pt, im_pt, s=dr_pt*10, c='r');
#plt.axis([0, 150, 0, 50])
plt.xlabel('Re')
plt.ylabel('Im')

grid(True)

savefig('/home/timy/Desktop/fig/crf_map.png')
plt.show()
