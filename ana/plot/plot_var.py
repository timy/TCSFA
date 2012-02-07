import matplotlib.pyplot as plt
import numpy as np
# this subroutine is used to find the relation between different physical quantities


traj_dir = '../../dat/'
selc_file_name = '../dat/selc_  1.dat'
n_traj_plot = 20

# read the info file
info = np.loadtxt( traj_dir + "info.dat" )

# get the index of the line
i_line = info[:,1]

# map the weight to the line width
w = info[:,2]
w_max = np.log10( np.amax(w) )
w_min = np.log10( np.amin(w) )
lw = 6.0 * ( np.log10( w ) - w_min ) / ( w_max - w_min ) + 0.1

# read the raw data
data = np.loadtxt( selc_file_name )


px_0, pz_0, ts_re, ts_im, x0, z0, px, pz, L, M_re, M_im, n_near_core, types = []
for i in range(n_traj_plot):
    px_0.append( data[i_line[i], 0] )
    pz_0.append( data[i_line[i], 1] )
    ts_re.append( data[i_line[i], 2] )
    ts_im.append( data[i_line[i], 3] )
    x0.append( data[i_line[i], 4] )
    z0.append( data[i_line[i], 5] )
    px.append( data[i_line[i], 6] )
    pz.append( data[i_line[i], 7] )
    L.append( data[i_line[i], 8] )
    M_re.append( data[i_line[i], 9] )
    M_im.append( data[i_line[i], 10] )
    n_near_core.append( data[i_line[i], 11] )
    types.append( data[i_line[i], 12] )
