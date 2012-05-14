import matplotlib.pyplot as plt
import numpy as np
import os

traj_dir = '../dat/traj/'
# n_traj_plot = len([f for f in os.listdir(traj_dir) 
#                    if f.startswith('traj_re') ])

# read the info file
info = np.loadtxt( traj_dir + "info.dat" )
n_traj_plot = sum(1 for line in open( traj_dir  + "info.dat" ))
print 'No. of trajectories: ', n_traj_plot


# map the weight to the line width
w = info[:,1]
w_max = np.log10( np.amax(w) )
w_min = np.log10( np.amin(w) )
lw = 6.0 * ( np.log10( w ) - w_min ) / ( w_max - w_min ) + 0.1
imax = np.argmax(w)
print("imax = %d" % imax)


# map the type to color
i_type = info[:,2]
cl = []
for i in range(n_traj_plot):
    if i_type[i] == 1:
        cl.append('r')
    elif i_type[i] == 2:
        cl.append('b')
    elif i_type[i] == 3:
        cl.append('g')
    elif i_type[i] == 4:
        cl.append('k')


# start ploting
plt.figure( figsize=(8,8), frameon=False )
ax = plt.subplot( 111, frameon=True )
plt.hold( True )

for i in range(n_traj_plot):
    data = np.loadtxt( traj_dir + "traj_re_%5d.dat"%(i+1) )
    plt.plot(data[:,2], data[:,3], lw=lw[i], label="%d"%i, c=cl[i] )
    plt.plot(data[0,2], data[0,3], 'd', c='r')

plt.xlim([-500,200])
plt.ylim([-20,10])
plt.grid(True)
plt.savefig( '../fig/traj.png' )
plt.show()
