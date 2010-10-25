import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *

file_name = 'energy_spec_qtm_all'
data = np.loadtxt('../data/'+file_name+'.dat')
data_energy = data[:,0]
data_angle = data[:,1]
data_hi = data[:,2]

n_grid_x = 250
n_grid_y = 100
grid_lower_x = 0
grid_upper_x = 360
grid_lower_y = 0.0
grid_upper_y = 0.5 * 27.2

#arry_pz = data_pz.reshape(n_grid_y, n_grid_x)
arry_energy = data_energy.reshape(n_grid_y, n_grid_x)
#arry_px = data_px.reshape(n_grid_y, n_grid_x)
arry_angle = data_angle.reshape(n_grid_y, n_grid_x)
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
fig = plt.figure(figsize=(16, 6), frameon = False )
ax = fig.add_subplot(111, aspect='auto')
fig.patch.set_alpha(0.0)
plt.imshow(np.rot90(arry_hi), interpolation='bilinear', 
           #origin='lower', #cmap=cm.gray,
           extent=[grid_lower_x, grid_upper_x,
                   grid_lower_y, grid_upper_y],
          aspect='auto')


#ax = gca()

#for tick in ax.xaxis.get_major_ticks():
#    tick.label1.set_fontsize(18)
#for tick in ax.yaxis.get_major_ticks():
#    tick.label1.set_fontsize(18)
#for tick in ax.xaxis.get_major_ticks():
#    tick.label1On = True
#    tick.label2On = False


#xlim(-2, 2)
#ylim(-1, 1)
clim(-3.8, -3)
#clim(250, 500)

#pcolor(arry_pz, arry_px, arry_hi)

colorbar()
grid()



savefig('/home/timy/Desktop/fig/'+file_name+'.svg')
plt.show()
