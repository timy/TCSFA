import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *

file_name = '../data/energy_spec_qtm_all'
data = np.loadtxt(file_name+'.dat')
data_energy = data[:,0]
data_angle = data[:,1]
data_w = data[:,2]

n_grid_energy = 250
n_grid_angle = 100
grid_lower_energy = 0.0
grid_upper_energy = 0.5
grid_lower_angle =  0.0
grid_upper_angle =  np.pi

arry_energy = data_energy.reshape(n_grid_angle, n_grid_energy)
arry_angle = data_angle.reshape(n_grid_angle, n_grid_energy)
arry_w = data_w.reshape(n_grid_angle, n_grid_energy)
ratio = amax(arry_w)
print 'minimum: ', amin(arry_w)
print 'maximum: ', amax(arry_w)
fig = plt.figure(figsize=(8, 6), frameon = False )
ax = fig.add_subplot(111)                           # aspect='equal')
fig.patch.set_alpha(0.0)


energy = []                  # x-axis
w = []                       # y-axis

for index in range(n_grid_energy):
    energy.append( arry_energy[0,index] * 27.2 )
    w.append( sum( arry_w[10:19,index] ) )

#w = log10(w)
plt.plot(energy, w, 'b-o')

plt.hold(True)

grid()



savefig('/home/timy/Desktop/fig/'+file_name+'.png')
plt.show()
