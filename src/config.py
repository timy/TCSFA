import ConfigParser
import sys
import os
from math import pi

def def_key_val( key_str, key_type, key_val ):
    if key_type == 'i':
        result_str = '#define ' + key_str + ' %15i\n'%key_val
    elif key_type == 'f':
        result_str = '#define ' + key_str + (' %15.8e\n' % key_val).replace('e', 'd')
    elif key_type == 's':
        result_str = '#define ' + key_str + '\t' + '\'' + key_val + '\'' + '\n'
    else:
        print 'Unknown data type of the key_val @ define_key_value'
        sys.exit()
    return result_str

def generate_lun():
    lun = 200
    while True:
        lun += 1
        yield lun
get_lun = generate_lun().next

config = ConfigParser.RawConfigParser()
config.read('../parameter.cfg')

dir_inc = './include/'
dir_ana_dat = os.getcwd() + '/../ana/dat/'
dir_raw_dat = os.getcwd() + '/../dat/'

if not os.path.exists( dir_inc ):
    os.makedirs( dir_inc )
if not os.path.exists( dir_raw_dat ):
    os.makedirs( dir_raw_dat )


######################### GRID ##################################
n_px     = config.getint('GRID', 'n_px')
n_pz     = config.getint('GRID', 'n_pz')
px_begin = config.getfloat('GRID', 'px_begin')
px_end   = config.getfloat('GRID', 'px_end')
pz_begin = config.getfloat('GRID', 'pz_begin')
pz_end   = config.getfloat('GRID', 'pz_end')

n_p0     = n_px*n_pz

if n_px == 1:
    d_px = 0.0
else:
    d_px     = (px_end - px_begin) / (n_px - 1.0)

if n_pz == 1:
    d_pz = 0.0
else:
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


with open(dir_inc+'inc_grid_size.h', 'w') as f:
    for s in grid_size:
        f.write( s )

with open(dir_inc+'inc_grid_boundary.h', 'w') as f:
    for s in grid_boundary:
        f.write( s )

########################## FIELD #################################
# field parameters
E0    = config.getfloat( 'FIELD', 'E0' )
om    = config.getfloat( 'FIELD', 'om' )
nc    = config.getfloat( 'FIELD', 'nc' )
xi    = config.getfloat( 'FIELD', 'xi' )
ph    = config.getfloat( 'FIELD', 'ph' )
envelope = config.get( 'FIELD', 'envelope' )

field = []
field.append( def_key_val( 'E0', 'f', E0 ) )
field.append( def_key_val( 'OM', 'f', om ) )
field.append( def_key_val( 'NC', 'f', nc ) )
field.append( def_key_val( 'XI', 'f', xi ) )
field.append( def_key_val( 'PH', 'f', ph ) )
field.append( def_key_val( 'Tp', 'f', 2*pi*nc/om))

if envelope=='sin2':
    field.append("#define PULSE_A_Z\t\t\tpulse_A_z_sin2\n")
    field.append("#define PULSE_E_Z\t\t\tpulse_E_z_sin2\n")
    field.append("#define PULSE_ALPHA_Z\tpulse_alpha_z_sin2\n")
    field.append("#define PULSE_ALPHA2_Z\tpulse_alpha2_z_sin2\n")
    field.append("#define PULSE_A_X\t\t\tpulse_A_x_sin2\n")
    field.append("#define PULSE_E_X\t\t\tpulse_E_x_sin2\n")
    field.append("#define PULSE_ALPHA_X\tpulse_alpha_x_sin2\n")
    field.append("#define PULSE_ALPHA2_X\tpulse_alpha2_x_sin2\n")
elif envelope=='const':
    field.append("#define PULSE_A_Z\t\t\tpulse_A_z_const\n")
    field.append("#define PULSE_E_Z\t\t\tpulse_E_z_const\n")
    field.append("#define PULSE_ALPHA_Z\tpulse_alpha_z_const\n")
    field.append("#define PULSE_ALPHA2_Z\tpulse_alpha2_z_const\n")
    field.append("#define PULSE_A_X\t\t\tpulse_A_x_const\n")
    field.append("#define PULSE_E_X\t\t\tpulse_E_x_const\n")
    field.append("#define PULSE_ALPHA_X\tpulse_alpha_x_const\n")
    field.append("#define PULSE_ALPHA2_X\tpulse_alpha2_x_const\n")

with open( dir_inc+'inc_field.h', 'w' ) as f:
    for s in field:
        f.write( s )

# field functions
field = []
if envelope == 'sin2':
    field.append("#include \'../pulse/pulse_sin2.f90\'\n")
elif envelope == 'const':
    field.append("#include \'../pulse/pulse_const.f90\'\n")

field.append( def_key_val( 'FID_PULSE', 'i', get_lun() ) )
if not os.path.exists(dir_ana_dat+'pulse/'):
    os.makedirs(dir_ana_dat+'pulse/')
field.append( def_key_val( 'FNM_PULSE', 's', dir_ana_dat+"pulse/pulse.dat") )

with open( dir_inc+'inc_field_func.h', 'w' ) as f:
    for s in field:
        f.write( s )

########################## ATOM ##################################
Ip      = config.getfloat( 'ATOM', 'Ip' )
charge  = config.getfloat( 'ATOM', 'charge' )

atom = []
atom.append( def_key_val( 'IONIZATION_IP', 'f', Ip ) )
atom.append( def_key_val( 'ATOM_CHARGE_Z', 'f', charge ) )

with open( dir_inc+'inc_atom.h', 'w' ) as f:
    for s in atom:
        f.write( s )

########################## TS_GUESS ##############################
tp = 2.0 * pi * nc / om

ts_max = config.getint( 'TS_GUESS', 'ts_max' )
lms_nx = config.getint( 'TS_GUESS', 'grid_re' )
lms_ny = config.getint( 'TS_GUESS', 'grid_im' )

if config.getboolean( 'TS_GUESS', 'auto' ):
    ts_re_lower = config.getfloat( 'TS_GUESS', 're_lower' )
    ts_re_upper = config.getfloat( 'TS_GUESS', 're_upper' )
    ts_re_lower = 0.01 + ts_re_lower * tp
    ts_re_upper = -0.01 + ts_re_upper * tp
else:
    ts_re_lower = config.getfloat( 'TS_GUESS', 're_lower' )
    ts_re_upper = config.getfloat( 'TS_GUESS', 're_upper' )
    if ts_re_lower < 0.0 or ts_re_lower > tp:
        sys.exit( '[TS_GUESS]: ts_re_lower invalid!' )
    if ts_re_upper < 0.0 or ts_re_upper > tp:
        sys.exit( '[TS_GUESS]: ts_re_upper invalid!' )

ts_im_lower = config.getfloat( 'TS_GUESS', 'im_lower' )
ts_im_upper = config.getfloat( 'TS_GUESS', 'im_upper' )

ts_guess = []
ts_guess.append( def_key_val( 'LMS_MAX_COUNT', 'i', ts_max ) )
ts_guess.append( def_key_val( 'LMS_NX       ', 'i', lms_nx ) )
ts_guess.append( def_key_val( 'LMS_NY       ', 'i', lms_ny ) )
ts_guess.append( def_key_val( 'LMS_RE_LOWER ', 'f', ts_re_lower ) )
ts_guess.append( def_key_val( 'LMS_RE_UPPER ', 'f', ts_re_upper ) )
ts_guess.append( def_key_val( 'LMS_IM_LOWER ', 'f', ts_im_lower ) )
ts_guess.append( def_key_val( 'LMS_IM_UPPER ', 'f', ts_im_upper ) )

with open( dir_inc+'inc_ts_guess.h', 'w' ) as f:
    for s in ts_guess:
        f.write( s )

############################# RK4 #################################
rk4_nt   = config.getint( 'RK4', 'nt' )
rk4_ne   = config.getint( 'RK4', 'ne' )
rk4_nmax = config.getint( 'RK4', 'nmax' )
rk4_eps  = config.getfloat( 'RK4', 'eps' )
rk4_r_threshold = config.getfloat( 'RK4', 'r_threshold' )


rk4 = []
rk4.append( def_key_val( 'RK4_NT           ', 'i', rk4_nt ) )
rk4.append( def_key_val( 'RK4_NE           ', 'i', rk4_ne ) )
rk4.append( def_key_val( 'RK4_NMAX         ', 'i', rk4_nmax ) )
rk4.append( def_key_val( 'RK4_EPS          ', 'f', rk4_eps ) )
rk4.append( def_key_val( 'RK4_R_THRESHOLD  ', 'f', rk4_r_threshold ) )


with open( dir_inc+'inc_rk4.h', 'w' ) as f:
    for s in rk4:
        f.write( s )

###################### REPRESENTATION ########################
representation_opt = config.get( 'MISC', 'representation' )

misc = [] 
misc.append( def_key_val( 'REPRESENTATION_OPT', 's', representation_opt ) )

###################### ANALYSIS & PLOT #######################
if len( sys.argv ) > 1:
    if 'traj' in sys.argv:
        print 'traj mode is selected'
        misc.append( def_key_val( 'MISC_PLOT_TRAJ ', 'i', 1 ) )
        if not os.path.exists(dir_ana_dat+'traj/'):
            os.makedirs(dir_ana_dat+'traj/')
        # plot trajectory/ies
        plot_rk4 = []
        plot_rk4.append( def_key_val( 'RK4_PLOT_TRAJ', 's', '') )
        plot_rk4.append( def_key_val( 'RK4_TRAJ_QUEU_FILE_ID', 'i', get_lun() ) )
        plot_rk4.append( def_key_val( 'RK4_TRAJ_QUEU_FILE_NAME', 's', dir_ana_dat+'traj/temp') )
        plot_rk4.append( def_key_val( 'RK4_TRAJ_PLOT_FILE_ID', 'i', get_lun() ) )
        plot_rk4.append( def_key_val( 'RK4_TRAJ_PLOT_FILE_NAME', 's', dir_ana_dat+'traj/traj') )
        
        with open( dir_inc+'inc_plot_rk4.h', 'w' ) as f:
            for s in plot_rk4:
                f.write( s )

    if 'crf' in sys.argv:
        print "crf mode is selected"
        misc.append( def_key_val( 'MISC_PLOT_CRF ', 'i', 1 ) )
        # trace the Complex-Root-Finding (CRF)
        plot_crf_track = []
        plot_crf_track.append( def_key_val( 'CRF_PLOT_STEP   ', 's', '' ) )
        plot_crf_track.append( def_key_val( 'CRF_PLOT_FILE_ID', 'i', get_lun() ) )
        plot_crf_track.append( def_key_val( 'CRF_PLOT_FILE_NAME', 's', dir_ana_dat+'crf_track.dat' ) )

        with open( dir_inc+'inc_plot_crf_track.h', 'w' ) as f:
            for s in plot_crf_track:
                f.write( s )

        # contour in the complex plane of sub-barrier region
        plot_crf_map = []
        plot_crf_map.append( def_key_val( 'CRF_PLOT_FILE_ID ', 'i', get_lun() ) )
        plot_crf_map.append( def_key_val( 'CRF_PLOT_N_GRID_X', 'i', lms_nx ) )
        plot_crf_map.append( def_key_val( 'CRF_PLOT_N_GRID_Y', 'i', lms_ny ) )
        plot_crf_map.append( def_key_val( 'CRF_PLOT_MAP', 's', '' ) )
        plot_crf_map.append( def_key_val( 'CRF_PLOT_MAP_FIX', 's', '' ) )
        plot_crf_map.append( def_key_val( 'CRF_PLOT_MAP_FILE_NAME', 's', dir_ana_dat+'crf_map.dat') )

        with open( dir_inc+'inc_plot_crf_map.h', 'w' ) as f:
            for s in plot_crf_map:
                f.write( s )

    if 'im_traj' in sys.argv:
        print "im_traj mode is selected"
        misc.append( def_key_val( 'MISC_PLOT_SUB ', 'i', 1 ) )
        # trace the integrand of the action
        plot_im_integrand = []
        plot_im_integrand.append( def_key_val( 'IM_PLOT_INTEGRAND  ', 's', '') )
        plot_im_integrand.append( def_key_val( 'IM_PLOT_FILE_ID', 'i', get_lun() ) )
        plot_im_integrand.append( def_key_val( 'IM_PLOT_FILE_NAME', 's', dir_ana_dat+'im_integrand.dat') )
        plot_im_integrand.append( def_key_val( 'IM_PLOT_N_PTS', 'i', 200 ) )
        plot_im_integrand.append( def_key_val( 'IM_PLOT_OFFSET', 'f', 40.0) )
        
        with open( dir_inc+'inc_plot_im_integrand.h', 'w' ) as f:
            for s in plot_im_integrand:
                f.write( s )

    if 'print' in sys.argv:
        print_info = 1
        misc.append( def_key_val( 'MISC_PRINT', 'i', 3 ) )

with open( dir_inc+'inc_misc.h', 'w' ) as f:
    for s in misc:
        f.write( s )

###################### CONSTANT  ########################
const = [] 
const.append( def_key_val( 'PI', 'f', pi ) )

with open( dir_inc+'inc_const.h', 'w' ) as f:
    for s in const:
        f.write( s )
