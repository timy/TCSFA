import matplotlib.pyplot as plt
import numpy as np

traj_dir = '../../dat/'
n_traj_plot = 400


# read the info file
info = np.loadtxt( traj_dir + "info.dat" )

# map the weight to the line width
w = info[:,2]
w_max = np.log10( np.amax(w) )
w_min = np.log10( np.amin(w) )
lw = 6.0 * ( np.log10( w ) - w_min ) / ( w_max - w_min ) + 0.1
imax = np.argmax(w)
print("imax = %d" % imax)


# map the type to color
i_type = info[:,3]
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
    data = np.loadtxt( traj_dir + "re_traj%5d.dat"%i )
    plt.plot(data[:,2], data[:,3], lw=lw[i], label="%d"%i, c=cl[i] )

plt.show()

