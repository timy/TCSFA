import ConfigParser
from math import pi

def def_key_val( key_str, key_type, key_val ):
    if key_type == 'i':
        result_str = '#define ' + key_str + ' %15i\n'%key_val
    elif key_type == 'f':
        result_str = '#define ' + key_str + (' %15.8e\n' % key_val).replace('e', 'd')
    elif key_type == 's':
        result_str = '#define ' + key_str + '\t' + key_val + '\n'
    else:
        print 'Unknown type of the key_val @ define_key_value'
        sys.exit()
    return result_str

config = ConfigParser.RawConfigParser()
config.read('../../parameter.cfg')

######################### GRID ##################################
n_px     = config.getint('GRID', 'n_px')
n_pz     = config.getint('GRID', 'n_pz')
px_begin = config.getfloat('GRID', 'px_begin')
px_end   = config.getfloat('GRID', 'px_end')
pz_begin = config.getfloat('GRID', 'pz_begin')
pz_end   = config.getfloat('GRID', 'pz_end')

n_p0     = n_px*n_pz
d_px     = (px_end - px_begin) / (n_px - 1.0)
d_pz     = (pz_end - pz_begin) / (n_pz - 1.0)

grid_size = []
grid_boundary = []

grid_size.append( def_key_val( 'N_PX', 'i', n_px ) )
grid_size.append( def_key_val( 'N_PZ', 'i', n_pz ) )
grid_size.append( def_key_val( 'N_P0', 'i', n_p0 ) )

grid_boundary.append( def_key_val( 'PX_BEGIN', 'f', px_begin ) )
grid_boundary.append( def_key_val( 'PX_END  ', 'f', px_end ) )
grid_boundary.append( def_key_val( 'PZ_BEGIN', 'f', pz_begin ) )
grid_boundary.append( def_key_val( 'PZ_END  ', 'f', pz_end ) )
grid_boundary.append( def_key_val( 'D_PX    ', 'f', d_px ) )
grid_boundary.append( def_key_val( 'D_PZ    ', 'f', d_pz ) )


with open('include/inc_grid_size.h', 'w') as f:
    for s in grid_size:
        f.write( s )

with open('include/inc_grid_boundary.h', 'w') as f:
    for s in grid_boundary:
        f.write( s )

########################## FIELD #################################
E0    = config.getfloat('FIELD', 'E0')
om    = config.getfloat('FIELD', 'om')
nc    = config.getfloat('FIELD', 'nc')
xi    = config.getfloat('FIELD', 'xi')

field = []
field.append( def_key_val( 'E0', 'f', E0 ) )
field.append( def_key_val( 'OM', 'f', om ) )
field.append( def_key_val( 'NC', 'f', nc ) )
field.append( def_key_val( 'XI', 'f', xi ) )

with open('include/inc_field.h', 'w') as f:
    for s in field:
        f.write( s )

########################## ATOM ##################################
Ip      = config.getfloat('ATOM', 'Ip')
charge  = config.getfloat('ATOM', 'charge')

atom = []
atom.append( def_key_val( 'IONIZATION_IP', 'f', Ip) )
atom.append( def_key_val( 'ATOM_CHARGE_Z', 'f', charge ) )

with open('include/inc_atom.h', 'w') as f:
    for s in atom:
        f.write( s )

########################## TS_GUESS ##############################
tp = 2.0 * pi * nc / om

ts_max = config.getint('TS_GUESS', 'ts_max')
lms_nx = config.getint('TS_GUESS', 'grid_re')
lms_ny = config.getint('TS_GUESS', 'grid_im')

if config.getboolean('TS_GUESS', 'auto'):
    ts_re_lower = 0.01
    ts_re_upper = tp-0.01
else:
    ts_re_lower = config.getfloat('TS_GUESS', 're_lower')
    ts_re_upper = config.getfloat('TS_GUESS', 're_upper')
    if ts_re_lower < 0.0 or ts_re_lower > tp:
        sys.exit('[TS_GUESS]: ts_re_lower invalid!')
    if ts_re_upper < 0.0 or ts_re_upper > tp:
        sys.exit('[TS_GUESS]: ts_re_upper invalid!')

ts_im_lower = config.getfloat('TS_GUESS', 'im_lower')
ts_im_upper = config.getfloat('TS_GUESS', 'im_upper')

ts_guess = []
ts_guess.append( def_key_val( 'LMS_MAX_COUNT', 'i', ts_max ) )
ts_guess.append( def_key_val( 'LMS_NX       ', 'i', lms_nx ) )
ts_guess.append( def_key_val( 'LMS_NY       ', 'i', lms_ny ) )
ts_guess.append( def_key_val( 'LMS_RE_LOWER ', 'f', ts_re_lower ) )
ts_guess.append( def_key_val( 'LMS_RE_UPPER ', 'f', ts_re_upper ) )
ts_guess.append( def_key_val( 'LMS_IM_LOWER ', 'f', ts_im_lower ) )
ts_guess.append( def_key_val( 'LMS_IM_UPPER ', 'f', ts_im_upper ) )

with open('include/inc_ts_guess.h', 'w') as f:
    for s in ts_guess:
        f.write( s )

############################# RK4 #################################
rk4_nt = config.getint('RK4', 'nt')
rk4_ne = config.getint('RK4', 'ne')
rk4_nmax = config.getint('RK4', 'nmax')
rk4_eps = config.getfloat('RK4', 'eps')

rk4 = []
rk4.append( def_key_val( 'RK4_NT  ', 'i', rk4_nt ) )
rk4.append( def_key_val( 'RK4_NE  ', 'i', rk4_ne ) )
rk4.append( def_key_val( 'RK4_NMAX', 'i', rk4_nmax ) )
rk4.append( def_key_val( 'RK4_EPS ', 'f', rk4_eps ) )

with open('include/inc_rk4.h', 'w') as f:
    for s in rk4:
        f.write( s )

































###################### CRF PLOT #######################
plot = config.getint('MISC', 'plot')

plot_crf_track = []
plot_crf_map = []
plot_rk4 = []

if plot == 0:

    plot_crf_track.append('')

if plot == 1:

    plot_crf_track.append( def_key_val( 'CRF_PLOT_STEP   ', 's', '' ) )
    plot_crf_track.append( def_key_val( 'CRF_PLOT_FILE_ID', 'i', 102 ) )
    plot_crf_track.append( def_key_val( 'CRF_PLOT_FILE_NAME', 's', '\'ana/data/crf_track.dat\'' ) )

    plot_crf_map.append( def_key_val( 'CRF_PLOT_FILE_ID ', 'i', 101 ) )
    plot_crf_map.append( def_key_val( 'CRF_PLOT_N_GRID_X', 'i', lms_nx ) )
    plot_crf_map.append( def_key_val( 'CRF_PLOT_N_GRID_Y', 'i', lms_ny ) )
    plot_crf_map.append( def_key_val( 'CRF_PLOT_MAP', 's', '' ) )
    plot_crf_map.append( def_key_val( 'CRF_PLOT_MAP_FIX', 's', '' ) )
    plot_crf_map.append( def_key_val( 'CRF_PLOT_MAP_FILE_NAME', 's', '\'ana/data/crf_map.dat\'') )

#define RK4_FILE_ID 101
#define RK4_PLOT_FILE_NAME 'dat/traj.dat'
    
    plot_rk4.append( def_key_val( 'RK4_PLOT_TRAJ', 's', '') )
    plot_rk4.append( def_key_val( 'RK4_TRAJ_QUEU_FILE_ID', 'i', 103 ) )
    plot_rk4.append( def_key_val( 'RK4_TRAJ_QUEU_FILE_NAME', 's', '\'ana/data/temp\'') )
    plot_rk4.append( def_key_val( 'RK4_TRAJ_PLOT_FILE_ID', 'i', 104 ) )
    plot_rk4.append( def_key_val( 'RK4_TRAJ_PLOT_FILE_NAME', 's', '\'ana/data/traj.dat\'') )

with open('include/inc_plot_crf_track.h', 'w') as f:
    for s in plot_crf_track:
        f.write( s )

with open('include/inc_plot_crf_map.h', 'w') as f:
    for s in plot_crf_map:
        f.write( s )

with open('include/inc_plot_rk4.h', 'w') as f:
    for s in plot_rk4:
        f.write( s )

