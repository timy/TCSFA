import numpy as np
from enthought.mayavi.mlab import *

data = np.loadtxt("../dat/series.dat");
#data = loadtxt("../dat/select 1.dat");
n_traj = len(data)
#n_traj = 0000
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

s = np.ones(n_traj) * 0.01

points3d(data_pz, data_px, data_ts_re, s, colormap="copper", scale_factor=.25)

