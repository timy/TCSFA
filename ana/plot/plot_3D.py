'''


'''

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib
import mpl_toolkits.mplot3d.axes3d as p3

data = np.loadtxt("../dat/series.dat");
#data = loadtxt("../dat/select 1.dat");
n_traj = len(data)
#n_traj = 1000
#n_traj = 55987
data_px_0 = data[:n_traj,0]
data_pz_0 = data[:n_traj,1]
data_ts_re = data[:n_traj,2]
data_ts_im = data[:n_traj,3]
data_x_0 = data[:n_traj,4]
data_z_0 = data[:n_traj,5]
data_px = data[:n_traj,6]
data_pz = data[:n_traj,7]
data_M_re = data[:n_traj,8]
data_M_im = data[:n_traj,9]
data_n_pass_x = data[:n_traj,10]
data_n_pass_z = data[:n_traj,11]
data_err = data[:n_traj,12]
data_type = data[:n_traj,13]

print 'number of trajectories:', n_traj

w = []
matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20) 

for i_traj in range(n_traj):
    w.append(data_M_re[i_traj]**2 + data_M_im[i_traj]**2)

ratio = 400 / ( np.amax( np.log10(w) ) - np.amin( np.log10(w) ) )
w_size = np.ceil( ratio * ( np.log10(w) - np.amin( np.log10(w) ) + 1 ) )

fig = plt.figure()
ax = p3.Axes3D(fig)

#ax = fig.add_subplot(111, projection='3d')

ax.scatter3D( data_pz_0, data_px_0, data_pz, s=1, c=data_pz )#data_ts_re/275.578, s=1, c=data_ts_re )
ax.view_init(60, 315)
#scatter( sqrt(data_pz**2+data_px**2), data_ts_re, c=data_type )#, s=w_size, c=data_type )
#scatter( data_pz_0, log10(data_M_re**2+data_M_im**2), c=data_type ) # s=w_size)
#xlim([410,420])
#ylim([10.2,10.8])
#plt.title( 't0 v.s. p_parallel ' + str( amax(w)) )
#plt.xlabel('p0_z')
#plt.ylabel('p0_x')


#for tick in ax.yaxis.get_major_ticks():
#    tick.label1On = False
#    tick.label2On = False
#ax.set_xticklabels(["", "", "", "", "", ""])
#ax.set_yticks([-40, -20, 0, 20, 40, 60, 80])
#ax.set_yticklabels(['', '', '', '', '', '', ''])
#ax.set_w_zticks([412,413,414,415,416,417,418])
#ax.set_ytickLabels(['','','','','','',''])

#for tick in ax.get_xticklabels():
#    tick.set_fontsize(2)
#for tick in ax.get_yticklabels():
#    tick.set_fontsize(2)

#grid(True)

plt.show()
