import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *

#data = np.loadtxt('../dat/spec 3.dat')
file_name = 'spec_qtm_123'
data = np.loadtxt('../data/'+file_name+'.dat')
data_pz = data[:,0]
data_px = data[:,1]
data_hi = data[:,2]

print amin(data_hi)

n_grid_x = 800 # should agree with the n_grid in file 'find_root.f90'
n_grid_y = 100


#transfer to 2D matrix
arry_pz = data_pz.reshape(n_grid_y, n_grid_x)
arry_px = data_px.reshape(n_grid_y, n_grid_x)
arry_hi = data_hi.reshape(n_grid_y, n_grid_x)
ratio = amax(arry_hi)
print ratio
#arry_hi = arry_hi / ratio
arry_hi = log10( arry_hi )

#arry_hi = arry_hi / ratio
#print amax(arry_hi)



#levels = arange(-3, 2, 0.2)
#cset1 = contourf(arry_re, arry_im, arry_hi, levels,
#                        cmap=cm.get_cmap('jet', len(levels)-1)
fig = plt.figure(figsize=(16,6), frameon = False )
ax = fig.add_subplot(111)
plt.imshow(arry_hi, interpolation='bilinear', #cmap=cm.gray,
                origin='lower', extent=[-2,2,0,0.5])


#ax = gca()

#for tick in ax.xaxis.get_major_ticks():
#    tick.label1.set_fontsize(18)
#for tick in ax.yaxis.get_major_ticks():
#    tick.label1.set_fontsize(18)
for tick in ax.xaxis.get_major_ticks():
    tick.label1On = True
    tick.label2On = False


xlim(-2, 2)
clim(-5.4, -1.4)
#clim(50, 150)

#pcolor(arry_pz, arry_px, arry_hi)

#colorbar()
grid()



savefig('/home/timy/Desktop/fig/'+file_name+'.png')
plt.show()
