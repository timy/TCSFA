import matplotlib.pyplot as plt
import numpy as np

file_name = 'hist_p0z.dat'
file_dir = '../dat/'

data = np.loadtxt( file_dir + file_name )

plt.plot( data[:,0], np.log10(data[:,1]), '*' )
plt.grid(True)
plt.show()
