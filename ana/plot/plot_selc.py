import matplotlib.pyplot as plt
import numpy as np
import os

file_name = 'spec_qtm_all.dat'
file_dir = '../dat/'
file_selc = '../proc/selc/sample.f'
file_fig = '../fig/selc.png'
nz, nx = 1600, 400
grid_lower_x, grid_upper_x = 0.0, 1.0
grid_lower_z, grid_upper_z = -2.0, 2.0
orders = 4.0

data = np.loadtxt( file_dir + file_name )
data_pz = data[:,0]
data_px = data[:,1]
data_w = data[:,2]

pz = data_pz.reshape( nx, nz )
px = data_px.reshape( nx, nz )
w = data_w.reshape( nx, nz )
w_max = np.amax(w)
w_lower_limit = np.log10(w_max) - orders
i_too_small = np.where( w < 10 ** w_lower_limit )
w[i_too_small] = 10 ** w_lower_limit
w = np.log10(w)


fig = plt.figure(figsize=(16,6), frameon = False )
ax = fig.add_subplot(111, aspect='equal')
fig.patch.set_alpha(0.0)
spec = plt.imshow( w, interpolation='bilinear',
            origin='lower', cmap=plt.cm.jet,
            extent=[grid_lower_z, grid_upper_z, grid_lower_x, grid_upper_x] )
plt.grid()


class SpectraAnalyzer:
    def __init__(self, imageSpectra):
        self.imag = imageSpectra.figure
        self.ptr_z = []
        self.ptr_x = []
        self.n_ptr = 0
        self.cid_click = imageSpectra.figure.canvas.mpl_connect(
            'button_press_event', self.on_click )
        self.cid_key = imageSpectra.figure.canvas.mpl_connect(
            'key_press_event', self.on_key )
    def on_click(self, event):
        print( 'button = %d, x = %d, y = %d, xdata = %f, ydata = %f' % 
               (event.button, event.x, event.y, event.xdata, event.ydata) )
        self.ptr_z.append(event.xdata)
        self.ptr_x.append(event.ydata)
        self.n_ptr = self.n_ptr + 1
        self.update()
    def on_key(self, event):
        print 'key', event.key, event.xdata, event.ydata
        if event.key == 'g':
            print 'start generating data file...'
            src_selc = []
            for i in range(self.n_ptr):
                src_selc.append('z0(%d) = %f\nx0(%d) = %f\n\n' % 
                                ((i+1), self.ptr_z[i], (i+1), self.ptr_x[i]))
            dir_fig = os.path.dirname( file_fig )
            if not os.path.exists( dir_fig ):
                os.makedirs( dir_fig )
            plt.savefig( file_fig )
            with open(file_selc, 'w') as f:
                for s in src_selc:
                    f.write( s )
            print 'done with generation.'
    def update(self):
        for i in range(self.n_ptr):
            plt.plot( self.ptr_z[i], self.ptr_x[i], '+', ms=20, c='k' )
            plt.text( self.ptr_z[i]+0.01, self.ptr_x[i]+0.01, "%d"%(i+1), 
                      color='k', fontsize=12)
            self.imag.canvas.draw()
        
        
sa = SpectraAnalyzer( spec )
plt.show()
