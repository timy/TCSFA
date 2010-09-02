'''
This file is used to plot the trajectory for a given p0

'''

from matplotlib import *
from pylab import *
import numpy as np




t0 = []
z0 = []
x0 = []
vz0 = []
vx0 = []
Mp_re = []
Mp_im = []
w = []
i_type = []

figure(figsize=(8,8), frameon=False)
ax = subplot(111, frameon=True)
hold(True)
criteria = np.loadtxt("../data/filter.dat")
traj_M = np.loadtxt("../data/traj_m.dat")
#n_traj = 514
n_traj = len(criteria)
index = range(n_traj)
#index = range(3,16)


print 'number of trajectories:', n_traj

# read in the amplitude and get ionization probability
for i_traj in index:
    data = np.loadtxt('../data/traj_m.dat');
    Mp_re.append( data[i_traj,0] );
    Mp_im.append( data[i_traj,1] );
    i_type.append( data[i_traj,2] );
    w.append( data[i_traj,0]**2 + data[i_traj,1]**2 );
    



print i_type
# read in the trajectories
for i_traj in index:
#    file_name = "../dat/traj_%3d.dat" % i_traj
    file_name = "../data/traj_%3d.dat" % criteria[i_traj]
    
    data = np.loadtxt(file_name)
    t = data[:,0];
    x = data[:,1];
    vx = data[:,2];
    z = data[:,3];
    vz = data[:,4];
    t0.append(t[0]);
    z0.append(z[0]);
    x0.append(x[0]);
    vz0.append(vz[0]);
    vx0.append(vx[0]);

    if i_type[i_traj] == 1:
        cl = 'r'
    elif i_type[i_traj] == 2:
        cl = 'b'
    elif i_type[i_traj] == 3:
        cl = 'g'
    elif i_type[i_traj] == 4:
        cl = 'k'
    plot( z, x, color=cl );

# plot the starting point for each trajectory
ratio = 400 / ( amax(log10(w)) - amin(log10(w)) )
w_size = ceil( ratio * ( log10( w ) - amin( log10(w) ) + 1 ) );
#print w_size
#scatter( t0, x0, c=t0, s=w_size )

#ax.set_xticks([-250, -200, -150, -100, -50, 0, 50, 100, 150, 200])
#ax.set_xticklabels(['', '', '', '', '', '', '', '', '', ''])
#ax.set_yticks([-40, -20, 0, 20, 40, 60, 80])
ax.set_yticklabels(['', '', '', '', '', '', ''])
grid(True)

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)



#plt.title( 'maximum of ionization probability: ' + str( amax(w)) )
#plt.xlabel('z')
#plt.ylabel('x')


plt.savefig('/data/Document/Poster/EMMI_Workshop_2010/latex-poster/fig/traj.svg')
plt.show()

