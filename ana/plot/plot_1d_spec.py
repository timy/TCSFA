import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors, ticker
from pylab import *


nz = 800
n_cut = 20;

data = np.loadtxt('../dat/spec_qtm_all.dat')
data_pz = data[:,0]
data_px = data[:,1]
data_hi = data[:,2]

pz = data_pz[0:nz];

spec = []

for i in range(nz):
    spec.append( 0.0 );

for i in range(nz):
    for j in range(n_cut):
        spec[i] = spec[i] + data_hi[nz*j+i];

print size(spec), size(pz)
plot( pz, spec )

xlabel('pz');
ylabel('ionization probability (Arb.)');
show()
