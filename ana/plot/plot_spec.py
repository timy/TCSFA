import matplotlib.pyplot as plt
import numpy as np

file_name = 'spec_qtm_all.dat'
file_dir = '../dat/'
nz, nx = 1600, 400
grid_lower_x, grid_upper_x = 0.0, 1.0
grid_lower_z, grid_upper_z = -2.0, 2.0
orders = 6.0


data = np.loadtxt( file_dir + file_name )
data_pz = data[:,0]
data_px = data[:,1]
data_w = data[:,2]

pz = data_pz.reshape( nx, nz )
px = data_px.reshape( nx, nz )
w = data_w.reshape( nx, nz )
w_max = np.amax(w)
print "w_max=", w_max
w_upper_limit = 10**2
w_lower_limit = np.log10(w_upper_limit) - orders
i_too_small = np.where( w < 10 ** w_lower_limit )
i_too_big = np.where( w > 10 ** w_upper_limit )
w[i_too_small] = 10 ** w_lower_limit
w[i_too_big] = 10 ** w_upper_limit
w = np.log10(w)


fig = plt.figure(figsize=(16,6), frameon = False )
ax = fig.add_subplot(111, aspect='equal')
fig.patch.set_alpha(0.0)
plt.imshow( w, interpolation='bilinear',
            origin='lower', #cmap=cm.gray,
            extent=[grid_lower_z, grid_upper_z, grid_lower_x, grid_upper_x] )
plt.grid()
plt.show()
