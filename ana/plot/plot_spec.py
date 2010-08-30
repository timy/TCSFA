import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *

file_name = 'spec_qtm_all'
data = np.loadtxt('../data/'+file_name+'.dat')
data_pz = data[:,0]
data_px = data[:,1]
data_hi = data[:,2]

n_grid_x = 500
n_grid_y = 300
grid_lower_x = -2.0
grid_upper_x = 2.0
grid_lower_y = -1.2
grid_upper_y = 1.2

arry_pz = data_pz.reshape(n_grid_y, n_grid_x)
arry_px = data_px.reshape(n_grid_y, n_grid_x)
arry_hi = data_hi.reshape(n_grid_y, n_grid_x)
ratio = amax(arry_hi)
print 'minimum: ', amin(arry_hi)
print 'maximum: ', amax(arry_hi)
#arry_hi = arry_hi / ratio
arry_hi = log10( arry_hi )

#arry_hi = arry_hi / ratio
#print amax(arry_hi)



#levels = arange(-3, 2, 0.2)
#cset1 = contourf(arry_re, arry_im, arry_hi, levels,
#                        cmap=cm.get_cmap('jet', len(levels)-1)
fig = plt.figure(figsize=(15, 8), frameon = False )
ax = fig.add_subplot(111)

plt.imshow(arry_hi, interpolation='bilinear', 
           origin='lower', #cmap=cm.gray,
           extent=[grid_lower_x, grid_upper_x,
                   grid_lower_y, grid_upper_y] )


#ax = gca()

#for tick in ax.xaxis.get_major_ticks():
#    tick.label1.set_fontsize(18)
#for tick in ax.yaxis.get_major_ticks():
#    tick.label1.set_fontsize(18)
for tick in ax.xaxis.get_major_ticks():
    tick.label1On = True
    tick.label2On = False


xlim(-2, 2)
ylim(-1.2, 1.2)
clim(-10, -1)
#clim(250, 500)

#pcolor(arry_pz, arry_px, arry_hi)

colorbar()
grid()



savefig('/home/timy/Desktop/fig/'+file_name+'.png')
plt.show()
