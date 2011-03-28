# ----------------------------------------------------------------------------
def plot_dist( filename, normalized=False ):
    from pylab import plot, plt
    import numpy as np

    data = np.loadtxt( '../../dat/' + filename );

    for i in range(1, 5):
        data[:,i] = np.log10(data[:,i]);
    if normalized:
        max_val = np.abs( np.amax( data[:,1:5] ) );
        for i in range(1, 5):
            data[:,i] = data[:,i] / max_val;

    plt.hold(True);
    plot( data[:,0], data[:,1], 'ro', label='T1', markersize=8 );
    plot( data[:,0], data[:,2], 'bo', label='T2', markersize=8 );
    plot( data[:,0], data[:,3], 'ko', label='T3', markersize=8 );
    plot( data[:,0], data[:,4], 'mo', label='T4', markersize=8 );
    plt.hold(False);

# ----------------------------------------------------------------------------
def plot_pulse( filename="pulse.dat", option="t-Ez", normalized=False, norm=False, style='r.-' ):
    from pylab import plot
    import numpy as np

    data = np.loadtxt('../../dat/' + filename );
     
    if normalized:
        for i in range(1, 4):
            max_val = np.amax( np.abs(data[:,i]) );
            data[:,i] = data[:,i] / max_val;
    if norm:
        for i in range(1, 4):
            data[:,i] = np.abs(data[:,i]);

    t, Az, Ax, Ez, Ex = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4];

    if option == 'Ez-Ex':
        return plot( Ez, Ex, style );
    elif option == 't-Ez':
        return plot( t, Ez, style, label='electric field' );
    elif option == 't-Ex':
        return plot( t, Ex, style );
    elif option == 't-Az':
        return plot( t, Az, style );
    elif option == 't-Ax':
        return plot( t, Ax, style );
    else:
        print "Option not available!";



def plot_traj( filename='traj.dat', option='z-x'):
    from pylab import plot, scatter
    import numpy as np

    data = np.loadtxt( '../data/' + filename );
    t, x, vx, z, vz = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4];

    plot( z, x );
    scatter( z[0], x[0], c='r', s=40)
