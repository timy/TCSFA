from tcsfa_plot import plot_dist, plot_pulse

# --------------------------------------------------------
def plot_dist_t0():
    return plot_dist( "dist_t0.dat", normalized=True );

# --------------------------------------------------------    
def plot_dist_z0():
    return plot_dist( "dist_z0.dat" );

# --------------------------------------------------------    
def plot_dist_L():
    return plot_dist( "dist_L.dat" );


from pylab import plt
from matplotlib.lines import Line2D
import matplotlib


matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)
fig = plt.figure(figsize=(10,6), frameon=False)
fig.patch.set_alpha(0.0)
plt_dist = plot_dist_t0();
plt.hold(True);
plt_pulse = plot_pulse(normalized=True);
plt.grid(True);
plt.xlim(250, 600);
plt.ylim(-2, 1.2);
plt.xlabel('t0 (a.u.)', fontsize=24, weight='bold');
plt.ylabel('Probability (Log10)', fontsize=24, weight='bold')

plt.legend(loc='best');
plt.title("pz=[-0.1,0.0], px=[0.0,0.2]", fontsize=24, weight='bold')
plt.savefig('/home/timy/Desktop/data/Document/Poster+Talks/DPG_2011/fig/t0_dist_(-0.3-0.2_0.05-0.15).png')
plt.show();
