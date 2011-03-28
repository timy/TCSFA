from tcsfa_plot import plot_dist, plot_pulse
from pylab import plt

plot_pulse( normalized=True )
plt.hold(True)
plot_pulse( option='t-Az', style='b', normalized=True )
plt.show()
