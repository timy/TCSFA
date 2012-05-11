import matplotlib.pyplot as plt
import numpy as np
import os

file_spec = 'spec_qtm_all.dat'
dir_spec = '../dat/'
dir_dat_selc = '../dat/selc/'
file_src_selc = '../proc/selc/sample.f'
file_src_selc_number = '../proc/selc/sample_number.f'
file_dat_selc = '../proc/selc/sample.dat'
file_fig = '../fig/selc.png'
nz, nx = 1600, 400
grid_lower_x, grid_upper_x = 0.0, 1.0
grid_lower_z, grid_upper_z = -2.0, 2.0
orders = 4.0

dz = 2.0 * (grid_upper_z - grid_lower_z) / (nz - 1.0)
dx = 2.0 * (grid_upper_x - grid_lower_x) / (nx - 1.0)

fig_lower_z = 1.0 * grid_lower_z
fig_upper_z = 1.0 * grid_upper_z
fig_lower_x = 1.0 * grid_lower_x
fig_upper_x = 1.0 * grid_upper_x

figsize_z, figsize_x = 7, 6
ratio_smp_lbl_z = 0.002 / figsize_z
ratio_smp_lbl_x = ratio_smp_lbl_z * 4

params = {#'backend': 'ps',
          'axes.labelsize': 20,
          'text.fontsize': 20,
          'legend.fontsize': 20,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'text.usetex': True,
          'figure.figsize': [figsize_z, figsize_x]}

plt.rcParams.update(params)

data = np.loadtxt( dir_spec + file_spec )
data_pz = data[:,0]
data_px = data[:,1]
data_w = data[:,2]

pz = data_pz.reshape( nx, nz )
px = data_px.reshape( nx, nz )
w = data_w.reshape( nx, nz )
#w_max = np.amax(w)
w_max = 1e0
w_lower_limit = np.log10(w_max) - orders
i_too_small = np.where( w < 10 ** w_lower_limit )
w[i_too_small] = 10 ** w_lower_limit
w = np.log10(w)


fig = plt.figure(figsize=(figsize_z,figsize_x), frameon = False )
ax = fig.add_subplot(111, aspect='equal')
fig.patch.set_alpha(0.0)
spec = plt.imshow( w, interpolation='bilinear',
            origin='lower', cmap=plt.cm.jet,
            extent=[fig_lower_z, fig_upper_z, fig_lower_x, fig_upper_x] )
#plt.xlim([-1.0, -0.8])
#plt.ylim([0.05, 0.25])
plt.xlabel(r'$p_z$ (a.u.)')
plt.ylabel(r'$p_x$ (a.u.)')
plt.grid()


class SpectraAnalyzer:
    def __init__(self, imageSpectra):
        self.imag = imageSpectra.figure
        self.smp_ptr_z, self.smp_ptr_x = [], []
        self.cur_ptr_z, self.cur_ptr_x = 0.0, 0.0
        self.n_smp = 0
        self.cid_click = imageSpectra.figure.canvas.mpl_connect(
            'button_press_event', self.on_click )
        self.cid_key = imageSpectra.figure.canvas.mpl_connect(
            'key_press_event', self.on_key )
#        self.cid_move = imageSpectra.figure.canvas.mpl_connect(
#            'motion_notify_event', self.on_move )
        self.cursor, = plt.plot( self.cur_ptr_z, self.cur_ptr_x, 'x', ms=20, color='k' )
#    def on_move(self, event):
#        if not event.inaxes: return
#        self.cur_ptr_z, self.cur_ptr_x = event.xdata, event.ydata
    def on_click(self, event):
        if not event.key == 'control': return
        if not event.inaxes: return
        print( 'button = %d, x = %d, y = %d, xdata = %f, ydata = %f' % 
               (event.button, event.x, event.y, event.xdata, event.ydata) )
        self.cur_ptr_z, self.cur_ptr_x = event.xdata, event.ydata
        self.mark()
    def on_key(self, event):
        print 'key', event.key, event.xdata, event.ydata
        
        if event.key == 'g':
            print 'start generating data file...'
            dir_fig = os.path.dirname( file_fig )
            if not os.path.exists( dir_fig ):
                os.makedirs( dir_fig )
            if not os.path.exists( dir_dat_selc ):
                os.makedirs( dir_dat_selc )
            # generate figure
            plt.savefig( file_fig )
            # generate source file for  FORTRAN
            src_selc = []
            for i in range(self.n_smp):
                src_selc.append('z0(%d) = %f\nx0(%d) = %f\n\n' % 
                                ((i+1), self.smp_ptr_z[i], (i+1), self.smp_ptr_x[i]))
            with open(file_src_selc, 'w') as f:
                for s in src_selc:
                    f.write(s)
            # write the the number of sample points
            src_selc_number = 'integer, parameter:: n_sel = %d' % self.n_smp
            with open(file_src_selc_number, 'w') as f:
                f.write(src_selc_number)
            # generate data file
            self.save()
            print 'done with generation.'
        if event.key == 'l':
            print 'start loading data from file \'sample.dat\'...'
            self.load()
            self.update()
        if event.key == 't':
            if event.inaxes:
                print 'the marker is targeted to the cursor...'
                self.cur_ptr_z, self.cur_ptr_x = event.xdata, event.ydata
                self.update_cursor()
        if event.key == 'm':
            self.mark()
        if event.key == 'left':
            pos = self.cur_ptr_z - dz
            if pos < grid_lower_z: return
            self.cur_ptr_z = pos
            self.update_cursor()
        if event.key == 'right':
            pos = self.cur_ptr_z + dz
            if pos > grid_upper_z: return
            self.cur_ptr_z = pos
            self.update_cursor()
        if event.key == 'up':
            pos = self.cur_ptr_x + dx
            if pos > grid_upper_x: return
            self.cur_ptr_x = pos
            self.update_cursor()
        if event.key == 'down':
            pos = self.cur_ptr_x - dx
            if pos < grid_lower_x: return
            self.cur_ptr_x = pos
            self.update_cursor()
    def mark(self):
            print 'the point is marked'
            self.smp_ptr_z.append(self.cur_ptr_z)
            self.smp_ptr_x.append(self.cur_ptr_x)
            self.n_smp = self.n_smp + 1
            self.update()
    def update(self):
        self.update_sample()
        self.update_cursor()
        print 'update done.'
    def update_sample(self):
        # get the size of the figure (display)
        len_fig_z, len_fig_x = plt.gcf().get_size_inches()
        pos_lbl_z = len_fig_z * ratio_smp_lbl_z
        pos_lbl_x = len_fig_x * ratio_smp_lbl_x
        # get the boundary of displayed data
#         box_dat = plt.axis()
#         len_dat_z = box_dat[1] - box_dat[0]
#         len_dat_x = box_dat[3] - box_dat[2]
#         pos_lbl_z = len_dat_z * ratio_smp_lbl
#         pos_lbl_x = len_dat_x * ratio_smp_lbl

        for i in range(self.n_smp):
            plt.plot( self.smp_ptr_z[i], self.smp_ptr_x[i], '+', ms=20, c='k' )
            plt.text( self.smp_ptr_z[i]+pos_lbl_z, self.smp_ptr_x[i]+pos_lbl_x, "%d"%(i+1), 
                      color='k', fontsize=12)
    def update_cursor(self):
        self.cursor.set_xdata([self.cur_ptr_z])
        self.cursor.set_ydata([self.cur_ptr_x])
        self.imag.canvas.draw()
    def load(self):
        ptr = np.loadtxt(file_dat_selc)
        self.n_smp = (ptr.shape)[0]
        self.smp_ptr_z = []
        self.smp_ptr_x = []
        for i in range(self.n_smp):
            self.smp_ptr_z.append(ptr[i,0])
            self.smp_ptr_x.append(ptr[i,1])
    def save(self):
        str_ptr = []
        for i in range(self.n_smp):
            str_ptr.append('%f %f\n' % (self.smp_ptr_z[i], self.smp_ptr_x[i]))
        with open(file_dat_selc, 'w') as f:
            for s in str_ptr:
                f.write(s)
        

sa = SpectraAnalyzer( spec )
plt.show()
