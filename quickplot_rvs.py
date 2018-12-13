#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Plot NRES RV results (distilled from PDF + table) for analysis.
#
# Rob Siverd
# Created:       2018-03-25
# Last modified: 2018-12-13
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.3.6"

## Python version-agnostic module reloading:
try:
    reload                              # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

## Modules:
import getopt
#import signal
#import glob
#import gc
import os
import sys
import time
#import ephem
import numpy as np
from numpy.lib.recfunctions import append_fields
#import datetime as dt
#import scipy.linalg as sla
#import scipy.signal as ssig
#import scipy.ndimage as ndi
#import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.ticker as mt
#import matplotlib._pylab_helpers as hlp
#from matplotlib.colors import LogNorm
#from matplotlib import colors
#import matplotlib.colors as mplcolors
import matplotlib.gridspec as gridspec
from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#np.set_printoptions(suppress=True, linewidth=160)
import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import PIL.Image as pli
#import seaborn as sns
import cmocean
import theil_sen as ts
#import window_filter as wf
#import itertools as itt

import lc
reload(lc)
kk = lc.LC()
kk.closeLCwin()
kk.lcwN = 5

import orbit
reload(orbit)

import qfit_orbit
reload(qfit_orbit)
qfo = qfit_orbit.OrbFit()

##--------------------------------------------------------------------------##

## Home-brew KDE:
#try:
#    import my_kde
#    reload(my_kde)
#    mk = my_kde
#except ImportError:
#    sys.stderr.write("\nError!  my_kde module not found!\n"
#           "Please install and try again ...\n\n")
#    sys.exit(1)

## Time conversion:
try:
    import astropy.time as astt
except ImportError:
    sys.stderr.write("\nError: astropy module not found!\n"
           "Please install and try again.\n\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
## Colors for fancy terminal output:
NRED    = '\033[0;31m'   ;  BRED    = '\033[1;31m'
NGREEN  = '\033[0;32m'   ;  BGREEN  = '\033[1;32m'
NYELLOW = '\033[0;33m'   ;  BYELLOW = '\033[1;33m'
NBLUE   = '\033[0;34m'   ;  BBLUE   = '\033[1;34m'
NMAG    = '\033[0;35m'   ;  BMAG    = '\033[1;35m'
NCYAN   = '\033[0;36m'   ;  BCYAN   = '\033[1;36m'
NWHITE  = '\033[0;37m'   ;  BWHITE  = '\033[1;37m'
ENDC    = '\033[0m'

## Suppress colors in cron jobs:
if (os.getenv('FUNCDEF') == '--nocolors'):
    NRED    = ''   ;  BRED    = ''
    NGREEN  = ''   ;  BGREEN  = ''
    NYELLOW = ''   ;  BYELLOW = ''
    NBLUE   = ''   ;  BBLUE   = ''
    NMAG    = ''   ;  BMAG    = ''
    NCYAN   = ''   ;  BCYAN   = ''
    NWHITE  = ''   ;  BWHITE  = ''
    ENDC    = ''

## Fancy text:
degree_sign = u'\N{DEGREE SIGN}'

## Dividers:
halfdiv = "----------------------------------------"
fulldiv = halfdiv + halfdiv

##--------------------------------------------------------------------------##
## Catch interruption cleanly:
#def signal_handler(signum, frame):
#    sys.stderr.write("\nInterrupted!\n\n")
#    sys.exit(1)
#
#signal.signal(signal.SIGINT, signal_handler)

##--------------------------------------------------------------------------##
def ldmap(things):
    return dict(zip(things, range(len(things))))

def argnear(vec, val):
    return (np.abs(vec - val)).argmin()

## Robust location/scale estimate using median/MAD:
def calc_ls_med_MAD(a, axis=None):
    """Return median and median absolute deviation of *a* (scaled to normal)."""
    med_val = np.median(a, axis=axis)
    sig_hat = (1.482602218 * np.median(np.abs(a - med_val), axis=axis))
    return (med_val, sig_hat)

## Robust location/scale estimate using median/IQR:
def calc_ls_med_IQR(a, axis=None):
    """Return median and inter-quartile range of *a* (scaled to normal)."""
    pctiles = np.percentile(a, [25, 50, 75], axis=axis)
    med_val = pctiles[1]
    sig_hat = (0.741301109 * (pctiles[2] - pctiles[0]))
    return (med_val, sig_hat)

## Select inliners given specified sigma threshold:
def pick_inliers(data, sig_thresh):
    med, sig = calc_ls_med_IQR(data)
    return ((np.abs(data - med) / sig) <= sig_thresh)

## Settings:
debug = False
timer = False
vlevel = 0
prog_name = 'quickplot_rvs.py'
full_prog = sys.argv[0]
base_prog = os.path.basename(full_prog)
num_todo = 0

## Options:
save_plot = None
ephem_file = None
min_npoints = 3
#zshift_file = None          # ZERO-specific velocity offsets
zero_offsets = {}           # ZERO file-specific velocity offsets
orbital_period = None

## Timekeeping:
day_zero = astt.Time("2018-01-01T00:00:00", scale='utc', format='isot')
#day_zero = 2458119.5    # 2018-01-01T00:00:00
#day_zero_string = "Jan 1, 2018"

##--------------------------------------------------------------------------##
## Argument type-checking:
def is_integer(asdf):
    try:
        int(asdf)
        return True
    except ValueError:
        return False

def is_float(asdf):
    try:
        float(asdf)
        return True
    except ValueError:
        return False

def parse_limit_string(lstring, ordered=True):
    """Returns lower,upper,valid tuple."""
    pieces = lstring.split(',')

    # Limit string must contain two elements:
    if len(pieces) != 2:
        return (None, None, False)      # exactly two elements required

    # Limits must be numbers:
    if not all([is_float(x) for x in pieces]):
        return (None, None, False)      # limits must be numbers

    # Extract:
    lower, upper = [float(x) for x in pieces]
    if ordered:
        lower, upper = (lower, upper) if (lower < upper) else (upper, lower)
    return lower, upper, True

##--------------------------------------------------------------------------##
##*********************     Help and options menu:     *********************##
##--------------------------------------------------------------------------##

## Syntax / how to run:
def usage(stream):
    stream.write("\n"
        + "Usage: %s [options] NRES_RV_data.txt\n" % base_prog
        + "Plot NRES RV results for inspection.\n"
        + "Version: %s\n" % __version__
        + "\n"
        + "Target information:\n"
        + "   -E, --ephems=FILE    data file with target ephemerides\n"
        + "   -P, --period=DAYS    target orbital period [days]\n"
        + "\n"
        + "Pipeline de-quirkification:\n"
        + "   -Z, --zshifts=FILE   load velocity shifts for ZERO files\n"
        + "\n"
        + "Available options:\n"
        + "       --debug          extra debugging info\n"
        + "   -h, --help           print this page\n"
        + "   -o, --output=FILE    save results to FILE\n"
        + "   -q, --quiet          suppress unnecessary output\n"
        + "   -t, --timer          report program run-time\n"
        + "   -v, --verbose        more status updates\n"
        + "\n")
        #+ "   -n, --numtodo=N     stop after N iterations\n"
        #+ "   -s, --sigcut=N      clip data beyond N*sigma\n"

##--------------------------------------------------------------------------##
##*********************       Parse command line:      *********************##
##--------------------------------------------------------------------------##

## Options:
short_opts = 'E:P:Z:o:hqtv' # n:s:
long_opts = ['ephems=', 'period=', 'zshifts=', 'output=',
                'debug', 'help', 'quiet', 'timer', 'verbose']
# 'numtodo=', 'sigcut='

## GNU-style parsing (with exception handling):
try:
    options, remainder = getopt.gnu_getopt(sys.argv[1:], short_opts, long_opts)
except getopt.GetoptError, err:
    sys.stderr.write("%s\n" % str(err))
    usage(sys.stderr)
    sys.exit(2)

## Handle selected options:
for opt, arg in options:
    # ------------------------------------------------
    if (opt == '--debug'):
        debug = True
        sys.stderr.write(BRED + "\nDebugging output enabled!" + ENDC + "\n")
    # ------------------------------------------------
    # ------------------------------------------------
    elif ((opt == '-E') or (opt == '--ephems')):
        if not os.path.isfile(arg):
            sys.stderr.write("File not found: %s\n" % arg)
            sys.exit(1)
        ephem_file = arg
        if (vlevel >= 0):
            msg = "Loading ephemerides from: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-P') or (opt == '--period')):
        if not is_float(arg):
            sys.stderr.write("Error!  Non-numeric orbit period: %s\n\n" % arg)
            sys.exit(1)
        orbital_period = float(arg)
        if (orbital_period <= 0):
            sys.stderr.write("Error: orbital period must be positive!\n")
            sys.exit(1)
        if (vlevel >= 0):
            msg = "Target orbital period: %.5f days" % orbital_period
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-Z') or (opt == '--zshifts')):
        if not os.path.isfile(arg):
            sys.stderr.write("Error, file not found: %s\n" % arg)
            sys.exit(1)
        vdata = np.genfromtxt(arg, dtype=None, names=True)
        zero_offsets.update(dict(zip(vdata['zfile'], vdata['vshift'])))
        #sys.stderr.write("vdata: %s\n" % str(vdata))
        if (vlevel >= 0):
            msg = "Loaded %d velocity shifts from: %s" % (len(vdata), arg)
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    #elif ((opt == '-n') or (opt == '--numtodo')):
    #    if not is_integer(arg):
    #        sys.stderr.write("Error!  Non-integer argument: %s\n\n" % arg)
    #        sys.exit(1)
    #    num_todo = int(arg)
    #    if (vlevel >= 0):
    #        msg = "Stopping after %d items." % num_todo
    #        sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-o') or (opt == '--output')):
        save_plot = arg
        if (vlevel >= 0):
            msg = "Saving RV plot to: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    #elif ((opt == '-s') or (opt == '--sigcut')):
    #    if not is_float(arg):
    #        sys.stderr.write("Error!  Non-numeric argument: %s\n\n" % arg)
    #        sys.exit(1)
    #    sigcut = float(arg)
    #    if (vlevel >= 0):
    #        msg = "Using %.2f-sigma outlier threshold." % sigcut
    #        sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-h') or (opt == '--help')):
        usage(sys.stdout)
        sys.exit(0)
    elif ((opt == '-q') or (opt == '--quiet')):
        vlevel -= 1
    elif ((opt == '-t') or (opt == '--timer')):
        timer = True
    elif ((opt == '-v') | (opt == '--verbose')):
        vlevel += 1
        sys.stderr.write(NYELLOW + "Increasing verbosity." + ENDC + "\n")
    # ------------------------------------------------
    else:
        msg = "Unhandled option: %s" % opt
        sys.stderr.write(BRED + "\n" + msg + ENDC + "\n\n")
        sys.exit(1)
    pass

## Verbosity:
if (vlevel >= 1):
    sys.stderr.write("%sVerbosity level: %d%s\n" % (NYELLOW, vlevel, ENDC))

## Full command line if highly verbose:
if (vlevel >= 2):
    sys.stderr.write("%s\nFull command line:%s\n" % (NCYAN, ENDC))
    sys.stderr.write("   %s\n" % sys.argv)

##--------------------------------------------------------------------------##
## Requirements check:
#requirements = [(bias_file, 'bias image'), (imlist_file, 'image list'),
#                (save_file, 'output file'), (fiber_pos[0], 'fiber X position'),
#                (fiber_pos[1], 'fiber Y position'),]
#for var,info in requirements:
#    if var == None:
#        sys.stderr.write(BRED + "No %s provided!" % info + ENDC + "\n")
#        usage(sys.stderr)
#        sys.exit(1)

## Output file is required:
if not save_plot:
    sys.stderr.write(BRED + "\nOutput file required!" + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Input file is required:
if (len(remainder) < 1):
    sys.stderr.write(BRED + "\nInput file required!" + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Check for input file:
data_file = remainder[0]
if not os.path.isfile(data_file):
    msg = "%s error:  file not found: %s" % (prog_name, data_file)
    sys.stderr.write("\n" + BRED + msg + ENDC + "\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Quick ASCII I/O:
#data_file = 'data.txt'
#all_data = np.genfromtxt(data_file, dtype=None)
raw_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True,
        delimiter=',')
raw_data = np.atleast_1d(raw_data)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True,
#                 delimiter='|', comments='%0%0%0%0')
#                 loose=True, invalid_raise=False)
#all_data = aia.read(data_file)
#all_data = pd.read_csv(data_file)
#all_data = pd.read_table(data_file, delim_whitespace=True)
#all_data = pd.read_table(data_file, sep='|')
fields = raw_data.dtype.names
n_raw_data = raw_data.size

## Numerical mapping for zero files:
unique_zero_files = list(set([x.rstrip('.fits') for x in raw_data['zero']]))
zero_file_mapping = ldmap(unique_zero_files)

## Apply offsets for known ZERO file offenders:
fix_cnames = ['rvval', 'ccval']
fix_data   = np.copy(raw_data)
use_offset = np.array([zero_offsets.get(x.rstrip('.fits'), 0.0) \
                        for x in fix_data['zero']])
for col in fix_cnames:
    fix_data[col] += use_offset

## Purge NaN-valued velocities:
vel_cnames = ['rvval', 'bcval', 'ccval']
velocities = np.column_stack([fix_data[x] for x in vel_cnames])
proc_fails = np.any(np.isnan(velocities), axis=1)
n_failures = np.sum(proc_fails)
all_data   = fix_data[(~proc_fails).nonzero()[0]]
n_all_data = all_data.size
sys.stderr.write("Valid data points: %d\n" % n_all_data)
use_data   = all_data

## Ignore data with low SNR:
#snr_cut = 10.0
#keepers = (all_data['snr'] >= snr_cut)
#all_data = all_data[keepers]

## Toss extreme outliers:
if (n_all_data >= 5):
    #keepers = pick_inliers(all_data['ccval'], 10)
    keepers = pick_inliers(all_data['rvval'], 10) \
            & (all_data['bjd_objs'] > 2450000.0)
    use_data = all_data[keepers]
else:
    sys.stderr.write("Skipping outlier rejection (too few points).\n")
n_use_data = len(use_data)
sys.stderr.write("Kept %d of %d data points.\n" % (n_use_data, n_raw_data))

## Stop if too few points to proceed:
if (n_use_data < min_npoints):
    sys.stderr.write("Insufficient data (need %d points)!\n" % min_npoints)
    sys.exit(1)

## Actionable quantities:
#use_tvec = timestamp - jd_offset
timestamp = use_data['bjd_obs']
jd_offset = timestamp.min()         # note for plotting
raw_velocity = use_data['rvval']

## Parse LCO site:
lco_site = np.array([x[:3] for x in use_data['tarball']])

## Numerical representation of zero file:
zf_number = np.array([zero_file_mapping[x] for x in use_data['zero']])

##--------------------------------------------------------------------------##
## Load ephemerides and look for match if file provided:
eph_data = []
targ_period = None
targ_eclmid = None
targ_t_peri = None
_have_ephem = False
nspec_phase = None
my_targname = None
if isinstance(ephem_file, str):
    eph_data = np.genfromtxt(ephem_file, dtype=None, 
            names=True, autostrip=True)
else:
    sys.stderr.write("No ephemeris file provided.\n")

## How to map ephemeris file to qfit_orbit parameters:
qfo_eph_cols = [('K', 'K'), ('Tp', 'tperi'), ('ecc', 'ecc'),
        ('per', 'period'), ('v0', 'v0'), ('w', 'omega')]

## Look for matching ephemeris:
targ_eph_params = {}
if (len(eph_data) > 0):
    for full_targ in use_data['target']:
        base_targ = full_targ.rstrip('_engr').rstrip('_bri')
        hits = (base_targ == eph_data['target'])
        if (np.sum(hits) == 1):
            eph_match = eph_data[hits]
            sys.stderr.write("Found ephemeris: %s\n" % str(eph_match))
            targ_t_peri = eph_match['tperi'][0]
            targ_period = eph_match['period'][0]
            targ_eclmid = eph_match['mid_eclipse'][0]
            my_targname = full_targ
            _have_ephem = True
            eph_match = eph_data[hits][0]
            for kq,kf in qfo_eph_cols:
                fval = eph_match[kf]
                #sys.stderr.write("kq: %s, kf: %s, fval: %f ... " % (kq,kf, fval))
                if not np.isnan(fval):
                #    sys.stderr.write("score!")
                    targ_eph_params[kq] = fval
                #else:
                #    sys.stderr.write("BOGUS!")
                #sys.stderr.write("\n")
            break
        pass
    pass

##--------------------------------------------------------------------------##
## Compute orbital phase:
if _have_ephem:
    #nspec_cycle = (timestamp - targ_t_peri) / targ_period
    nspec_cycle = (timestamp - targ_eph_params['Tp']) / targ_eph_params['per']
    nspec_phase = nspec_cycle % 1.0
else:
    nspec_phase = np.zeros_like(timestamp)

##--------------------------------------------------------------------------##
## Timestamp modification:
def time_warp(jdutc, jd_offset, scale):
    return (jdutc - jd_offset) * scale

## Self-consistent time-modification for plotting:
#tfudge = partial(time_warp, jd_offset=jd_offset, scale=1.0)    # relative days
tfudge = partial(time_warp, jd_offset=day_zero.jd, scale=1.0)    # relative days
#tfudge = partial(time_warp, jd_offset=tstart.jd, scale=24.0)    # relative hrs
#tfudge = partial(time_warp, jd_offset=tstart.jd, scale=1440.0)  # relative min

##--------------------------------------------------------------------------##
## Misc:
#def log_10_product(x, pos):
#   """The two args are the value and tick position.
#   Label ticks with the product of the exponentiation."""
#   return '%.2f' % (x)  # floating-point
#
#formatter = plt.FuncFormatter(log_10_product) # wrap function for use

## Convenient, percentile-based plot limits:
#def nice_limits(vec, pctiles=[1,99], pad=1.2):
#    ends = np.percentile(vec, pctiles)
#    middle = np.average(ends)
#    return (middle + pad * (ends - middle))

## Convenient plot limits for datetime/astropy.Time content:
#def nice_time_limits(tvec, buffer=0.05):
#    lower = tvec.min()
#    upper = tvec.max()
#    ndays = upper - lower
#    return ((lower - 0.05*ndays).datetime, (upper + 0.05*ndays).datetime)

## Convenient limits for datetime objects:
#def dt_limits(vec, pad=0.1):
#    tstart, tstop = vec.min(), vec.max()
#    trange = (tstop - tstart).total_seconds()
#    tpad = dt.timedelta(seconds=pad*trange)
#    return (tstart - tpad, tstop + tpad)

##--------------------------------------------------------------------------##
## Velocity offsets (date, shift):
#vel_offset_spec = [(0.0, 0.0), (2458119.5, 0.0), (2458219.5, 60.0)]
#
#zero_point_offsets = []
#for rvjd,rvvel in zip(timestamp, raw_velocity):
#    offset = [vshift for vsjd,vshift in vel_offset_spec if rvjd >= vsjd][-1]
#    zero_point_offsets.append(offset)
#zero_point_offsets = np.array(zero_point_offsets)
#
#nspec_radvel = raw_velocity + zero_point_offsets
nspec_radvel = raw_velocity

##--------------------------------------------------------------------------##
## Theil-Sen line-fitting for linear trend removal:
#icept, slope = ts.linefit(timestamp, nspec_radvel)
trend = ts.linefit(timestamp, nspec_radvel)
sys.stderr.write("Got linear fit of %10.5f km/s/day\n" % trend[1])

trend[1] += 0.5
#detrend_RV = nspec_radvel - (trend[0] + timestamp * trend[1])
detrend_RV = nspec_radvel - (timestamp * trend[1])
detrend_RV -= np.median(detrend_RV)


##--------------------------------------------------------------------------##
##------------------         Pick a Velocity to Use         ----------------##
##--------------------------------------------------------------------------##

use_detrended = False
#use_detrended = True
#if use_detrended:
#    nspec_radvel -= trend[0] + timestamp * trend
star_vel = detrend_RV if use_detrended else nspec_radvel

##--------------------------------------------------------------------------##
##------------------           Basic RV Analysis            ----------------##
##--------------------------------------------------------------------------##

if _have_ephem:
    #qfo.set_data(timestamp, nspec_radvel)
    qfo.set_data(timestamp, star_vel)
    #qfo.set_element('per', targ_period) # known orbital period
    #qfo.set_element('Tp', targ_t_peri)  # guess Tperi ~= Tc
    ##qfo.set_element('Tp', targ_t_peri-1.0)
    ##qfo.set_element('Tp', targ_t_peri + 0.5*targ_period)  # guess Tperi ~= Tc
    #qfo.set_element('ecc', 0.01)
    #pshift = 0.0

    # Attempt blind (awesome) orbital fit:
    #sfargs = {eccmax=0.95}
    #sfargs = {'n_ecc':30, 'n_omega':30, 'maxdepth':2}
    sfargs = {'n_ecc':20, 'n_omega':20, 'maxdepth':1, 'e_max':0.95}
    #sfargs = {'n_ecc':20, 'n_omega':20, 'maxdepth':1, 'full':True}
    #sfargs = {'n_ecc':20, 'n_omega':20, 'maxdepth':2}
    #(ecc, ww, pshift, resid), gdata = qfo.shapefit(nspec_phase, nspec_radvel, **sfargs)
    #ecc, ww, pshift, resid = qfo.shapefit(nspec_phase, nspec_radvel, **sfargs)
    ecc, ww, pshift, resid = qfo.shapefit(nspec_phase, star_vel, **sfargs)
    #ecc, ww, pshift, resid = qfo.shapefit(nspec_phase, nspec_radvel)
    qfo.set_elements({'ecc':ecc, 'w':ww})

    # MANUAL OVERRIDE:
    ##qfo.set_elements({'ecc':0.308, 'w':np.radians(166.), 'Tp':2455583.722984})
    #qfo.set_elements({'ecc':0.320})
    ##qfo.set_elements({'w':np.radians(166.)})
    #qfo.set_elements({'w':np.radians(220.)})
    ##qfo.set_elements({'Tp':2455583.800})
    #qfo.set_elements({'K':39.0})
    #qfo.set_elements({'v0': -102.0})
    qfo.set_elements(targ_eph_params)
    #pshift = 0.55
    orb_par_txt = '\n'.join(["%3s: %.4f" % P for P in qfo.elem_dict.items()])

    #if (my_targname == 'KEBC09C05946'):
    #    qfo.set_element('ecc', 0.5)
    #    qfo.set_element('w', 4.10)
    #    #pshift = 0.55
    #    #pshift = 0.0
    #    pshift = 0.45

#omegas = 2.0 * np.pi * np.linspace(0.0, 1.0, 10, endpoint=False, dtype='float')
#eccval = 0.8
#for ww in omegas:
#    phi, vel = orbit.orbit_curve(eccval, ww, norm=True)
#    #rvmin, rvmax = curve.min(), curve.max()
#    sys.stderr.write("For ww=%5.3f, RV min/max = %8.4f, %8.4f\n"
#            % (ww, vel.min(), vel.max()))


##--------------------------------------------------------------------------##
## Group data points by site:
tf = pd.DataFrame(dict(time=tfudge(timestamp), phase=nspec_phase,
                    rv=star_vel, site=lco_site, zfile=zf_number))
                    #rv=nspec_radvel, site=lco_site))
sitegroups = tf.groupby('site')

##--------------------------------------------------------------------------##
## Plot config:

# gridspec examples:
# https://matplotlib.org/users/gridspec.html

#gs1 = gridspec.GridSpec(4, 4)
#gs1 = gridspec.GridSpec(2, 2)
if _have_ephem:
    gs1 = gridspec.GridSpec(2, 1)
else:
    gs1 = gridspec.GridSpec(1, 1)
#gs1.update(wspace=0.025, hspace=0.05)  # set axis spacing

##--------------------------------------------------------------------------##
fig_dims = (12, 10)
fig = plt.figure(1, figsize=fig_dims)
plt.gcf().clf()
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = fig.add_subplot(211)
ax1 = plt.subplot(gs1[0, 0])
#ax1 = fig.add_axes([0, 0, 1, 1])
#ax1.patch.set_facecolor((0.8, 0.8, 0.8))
ax1.grid(True)
dotkw = {'marker':'o', 'ls':'', 'ms':4}
#ax1.set_xlabel
for site, group in sitegroups:
    #ax1.scatter(group.time, group.rv, s=50, c=group.zfile)
    ax1.plot(group.time, group.rv, label=site, **dotkw)
#ax1.set_xlabel("Days since %.3f" % jd_offset)
ax1.set_xlabel("Days since %.3f (%s)" % (day_zero.jd, day_zero.iso))
ax1.set_ylabel("Raw RV (km/s)")
ax1.legend(loc='best')

##ax1.axis('off')
##ax2 = fig.add_subplot(212) ax2 = plt.subplot(gs1[1, 0])
#ax2 = plt.subplot(gs1[1, 0])
#ax2.grid(True)
#ax2.scatter(tfudge(timestamp), detrend_RV, lw=0, s=25)
#ax2.set_xlabel("Days since %.3f" % jd_offset)
#ax2.set_ylabel("Detrended RV (km/s)")

if _have_ephem:
    #ax3 = plt.subplot(gs1[0, 1])
    ax3 = plt.subplot(gs1[1, 0])
    ax3.grid(True)
    for site, group in sitegroups:
        #ax3.scatter(group.phase, group.rv, s=50, c=group.zfile)
        ax3.plot(group.phase, group.rv, label=site, **dotkw)
    ax3.set_xlim(-0.05, 1.05)
    ax3.set_xlabel("Orbital Phase (P=%.5f days)" % targ_period)
    ax3.set_ylabel("Raw RV (km/s)")
    tphase, tcurve = qfo.make_orbit(points=2000, norm=True)
    #tphase, tcurve = qfo.test_orbit(points=1000, norm=True)
    aphase = (tphase - pshift) % 1.0
    #ax3.scatter(tphase, tcurve, lw=0, s=1, c='g', label='not shifted')
    #ax3.scatter(aphase, tcurve, lw=0, s=1, c='r', label='shifted %.3f'%pshift)
    ax3.scatter(aphase, tcurve, lw=0, s=1, c='r', label='shapefit')
    #ax3.plot(tphase, tcurve, lw=1, ls=':', c='r')

    orb_par_txt = '\n'.join(["%3s: %.4f" % P for P in qfo.elem_dict.items()])
    ax3.text(0.99, 0.01, orb_par_txt, transform=ax3.transAxes, ha='right',
            va='bottom', multialignment='right',
            fontdict={'family':'monospace'},
            bbox=dict(facecolor='white', pad=5.0))
    ax3.legend(loc='best')

    #shift,resid = qfo._test_orbit_shape(nspec_phase, nspec_radvel,
    shift,resid = qfo._test_orbit_shape(nspec_phase, star_vel,
            tphase, tcurve)

    # estimate scatter relative to curve:
    rv_diffs = []
    for nphi,nvel in zip(nspec_phase, star_vel):
        #sys.stderr.write("nphi: %f, nvel: %f\n" % (nphi, nvel))
        #which = argnear(aphase, nphi)
        nearest_vel = tcurve[argnear(aphase, nphi)]
        #sys.stderr.write("which: %d\n" % which)
        #sys.stderr.write("near_vel: %f\n" % nearest_vel)
        rv_diffs.append(nvel - nearest_vel)
    rv_diffs = np.array(rv_diffs)
    useful = pick_inliers(rv_diffs, 4.0)

#if _have_ephem:
#    ax4 = plt.subplot(gs1[1, 1])
#    ax4.grid(True)
#    ax4.scatter(nspec_phase, detrend_RV, lw=0, s=25)
#    ax4.set_xlim(-0.05, 1.05)
#    ax4.set_xlabel("Orbital Phase (P=%.5f days)" % targ_period)
#    ax4.set_ylabel("Detrended RV (km/s)")

## Disable axis offsets:
#ax1.xaxis.get_major_formatter().set_useOffset(False)
#ax1.yaxis.get_major_formatter().set_useOffset(False)

#ax1.plot(kde_pnts, kde_vals)

#blurb = "some text"
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes)
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes,
#      va='top', ha='left', bbox=dict(facecolor='white', pad=10.0))

#colors = cm.rainbow(np.linspace(0, 1, len(plot_list)))
#for camid, c in zip(plot_list, colors):
#    cam_data = subsets[camid]
#    xvalue = cam_data['CCDATEMP']
#    yvalue = cam_data['PIX_MED']
#    yvalue = cam_data['IMEAN']
#    ax1.scatter(xvalue, yvalue, color=c, lw=0, label=camid)

#mtickpos = [2,5,7]
#ndecades = 1.0   # for symlog, set width of linear portion in units of dex
#nonposx='mask' | nonposx='clip' | nonposy='mask' | nonposy='clip'
#ax1.set_xscale('log', basex=10, nonposx='mask', subsx=mtickpos)
#ax1.set_xscale('log', nonposx='clip', subsx=[3])
#ax1.set_yscale('symlog', basey=10, linthreshy=0.1, linscaley=ndecades)
#ax1.xaxis.set_major_formatter(formatter) # re-format x ticks
#ax1.set_ylim(ax1.get_ylim()[::-1])
#ax1.set_xlabel('whatever', labelpad=30)  # push X label down 

#ax1.set_xticks([1.0, 3.0, 10.0, 30.0, 100.0])
#ax1.set_xticks([1, 2, 3], ['Jan', 'Feb', 'Mar'])
#for label in ax1.get_xticklabels():
#    label.set_rotation(30)

#ax1.set_xlim(nice_limits(xvec, pctiles=[1,99], pad=1.2))
#ax1.set_ylim(nice_limits(yvec, pctiles=[1,99], pad=1.2))

#spts = ax1.scatter(x, y, lw=0, s=5)
#cbar = fig.colorbar(spts, orientation='vertical')
#cbar.formatter.set_useOffset(False)
#cbar.update_ticks()

fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
plt.draw()

if (save_plot != None):
    fig.savefig(save_plot)





######################################################################
# CHANGELOG (quickplot_rvs.py):
#---------------------------------------------------------------------
#
#  2018-06-07:
#     -- Increased __version__ to 0.3.5.
#     -- Added support for ZERO-specific RV offsets.
#
#  2018-06-06:
#     -- Increased __version__ to 0.3.0.
#     -- Now starting to keep track of systematic shifts in data.
#     -- Unphased RV plots (top panel) now use the same date as a reference
#           point so that data can be compared among objects. This should make
#           it much easier to spot systematic changes and issues.
#
#  2018-05-31:
#     -- Increased __version__ to 0.2.5.
#     -- Now support orbital parameter retrieval from ephems file.
#
#  2018-05-15:
#     -- Increased __version__ to 0.2.1.
#     -- Experiments with trend removal.
#
#  2018-03-27:
#     -- Increased __version__ to 0.2.0.
#     -- Now have automated shape-based orbit estimation!
#
#  2018-03-26:
#     -- Increased __version__ to 0.1.0.
#     -- Introduced orbit.py for simple fitting.
#
#  2018-03-25:
#     -- Increased __version__ to 0.0.5.
#     -- Basic plots working!
#     -- Increased __version__ to 0.0.1.
#     -- First created quickplot_rvs.py.
#
