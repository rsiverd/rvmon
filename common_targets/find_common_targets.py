#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Identify which objects have been observed by multiple NRES units.
#
# Rob Siverd
# Created:       2018-08-09
# Last modified: 2018-08-09
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.0.1"

## Python version-agnostic module reloading:
try:
    reload                              # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

## Modules:
#import argparse
import os
import sys
import time
import numpy as np
#from numpy.lib.recfunctions import append_fields
#import datetime as dt
#from dateutil import parser as dtp
#import scipy.linalg as sla
#import scipy.signal as ssig
#import scipy.ndimage as ndi
#import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd
#from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#np.set_printoptions(suppress=True, linewidth=160)
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import PIL.Image as pli
#import seaborn as sns
#import cmocean
#import theil_sen as ts
#import window_filter as wf
#import itertools as itt

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

## Fast FITS I/O:
#try:
#    import fitsio
#except ImportError:
#    sys.stderr.write("\nError: fitsio module not found!\n")
#    sys.exit(1)

## FITS I/O:
#try:
#    import astropy.io.fits as pf
#except ImportError:
#    try:
#       import pyfits as pf
#    except ImportError:
#        sys.stderr.write("\nError!  No FITS I/O module found!\n"
#               "Install either astropy.io.fits or pyfits and try again!\n\n")
#        sys.exit(1)

## ASCII I/O:
#try:
#    import astropy.io.ascii as aia
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n")
#    sys.exit(1)

## Time conversion:
#try:
#    import astropy.time as astt
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n"
#           "Please install and try again.\n\n")
#    sys.exit(1)

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




##--------------------------------------------------------------------------##
## Parse arguments and run script:
#class MyParser(argparse.ArgumentParser):
#    def error(self, message):
#        sys.stderr.write('error: %s\n' % message)
#        self.print_help()
#        sys.exit(2)
#
#if __name__ == '__main__':
#
#    # ------------------------------------------------------------------
#    descr_txt = """
#    PUT DESCRIPTION HERE.
#    
#    Version: %s
#    """ % __version__
#    parser = argparse.ArgumentParser(
#            prog='PROGRAM_NAME_HERE',
#            prog=os.path.basename(__file__),
#            description='PUT DESCRIPTION HERE.')
#            #description=descr_txt)
#    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt)
#    parser.add_argument('firstpos', help='first positional argument')
#    parser.add_argument('-s', '--site',
#            help='Site to retrieve data for', required=True)
#    parser.add_argument('-n', '--number_of_days', default=1,
#            help='Number of days of data to retrieve.')
#    parser.add_argument('-o', '--output_file', 
#            default='observations.csv', help='Output filename.')
#    parser.add_argument("--start", type=str, default=None, 
#            help="Start time for date range query.")
#    parser.add_argument("--end", type=str, default=None,
#            help="End time for date range query.")
#    parser.add_argument('-d', '--dayshift', required=False, default=0,
#            help='Switch between days (1=tom, 0=today, -1=yest', type=int)
#    parser.add_argument('-e', '--encl', nargs=1, required=False,
#            help='Encl to make URL for', choices=all_encls, default=all_encls)
#    parser.add_argument('-s', '--site', nargs=1, required=False,
#            help='Site to make URL for', choices=all_sites, default=all_sites)
#    parser.add_argument('-q', '--quiet', action='count', default=0,
#            help='less progress/status reporting')
#    parser.add_argument('-v', '--verbose', action='count', default=0,
#            help='more progress/status reporting')
#    parser.add_argument('--debug', dest='debug', default=False,
#            help='Enable extra debugging messages', action='store_true')
#    parser.add_argument('remainder', help='other stuff', nargs='*')
#    # ------------------------------------------------------------------
#    # ------------------------------------------------------------------
#    ofgroup = parser.add_argument_group('Output format')
#    fmtparse = ofgroup.add_mutually_exclusive_group()
#    fmtparse.add_argument('--python', required=False, dest='output_mode',
#            help='Return Python dictionary with results [default]',
#            default='pydict', action='store_const', const='pydict')
#    bash_var = 'ARRAY_NAME'
#    bash_msg = 'output Bash code snippet (use with eval) to declare '
#    bash_msg += 'an associative array %s containing results' % bash_var
#    fmtparse.add_argument('--bash', required=False, default=None,
#            help=bash_msg, dest='bash_array', metavar=bash_var)
#    # ------------------------------------------------------------------
#
#    context = parser.parse_args()
#    #context.vlevel = context.verbose - context.quiet
#    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)
#
##--------------------------------------------------------------------------##
## Object name extractor:
def get_clean_objname(objstring):
    suffixes = ['_pystrat', '_engr', '_bri']
    fib0, fib2 = objstring.split('&thar&')
    fib_names = [x for x in objstring.split('&thar&') if (x != 'none')]
    if (len(fib_names) != 1):
        return "parse_error"
    clean_name = fib_names[0]
    for suff in suffixes:
        clean_name = clean_name.rstrip(suff)
    clean_name = clean_name.replace(' ', '')
    #for i in range(len(clean_name)):
    #    clean_name = clean_
    return clean_name

def get_site(filename):
    return filename.split('/')[3]

##--------------------------------------------------------------------------##
## Quick ASCII I/O:
data_file = 'observed.txt'
#all_data = np.genfromtxt(data_file, dtype=None)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True)
all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True,
                 delimiter=',', comments='%0%0%0%0')
#                 loose=True, invalid_raise=False)
#all_data = aia.read(data_file)

## Make sure objects string has only 2 ampersands:
#bogus = np.bool_([(len(x.split('&')) != 3) for x in all_data['objects']])
bogus = np.bool_([(len(x.split('&thar&')) != 2) for x in all_data['objects']])
use_data = all_data[~bogus]

#for ii in range(10, 20):
#    dirty_name = use_data['objects'][ii]
#    clean_name = get_clean_objname(dirty_name)
#    sys.stderr.write("%s --> %s\n" % (dirty_name, clean_name))

## Sanitize object names and identify unique targets:
clean_site_list = np.array([get_site(x) for x in use_data['filename']])
clean_obj_names = np.array([get_clean_objname(x) for x in use_data['objects']])
known_objects = np.unique(clean_obj_names)

## Count number of files per clean object name:
target_imcount = \
        dict([(x, np.sum(clean_obj_names == x)) for x in known_objects])

## Identify objects with multiple spectra:
repeated_targs = np.array([kk for kk,vv in target_imcount.items() if (vv > 1)])

## Count permutations of site + OBJECTS for each sanitized target name:
for target in repeated_targs:
    which = (clean_obj_names == target)
    total = np.sum(which)
    #sys.stderr.write("%s --> %d hits\n" % (target, total))
    subset = use_data[which]
    #sites = clean_site_list[which]
    descr_list = []
    for stuff in subset:
        lsite = get_site(stuff['filename'])
        sys.stderr.write("%s %s\n" % (lsite, stuff['objects']))
        descr_list.append("%s_%s" % (lsite, stuff['objects']))
    permutations = set(descr_list)
    #sys.stderr.write("permutations: %s\n" % str(descr_list))
    if (len(permutations) < 2):
        continue    # not useful
    sys.stderr.write("n_permutations: %d\n" % len(permutations))
    

##--------------------------------------------------------------------------##


######################################################################
# CHANGELOG (find_common_targets.py):
#---------------------------------------------------------------------
#
#  2018-08-09:
#     -- Increased __version__ to 0.0.1.
#     -- First created find_common_targets.py.
#
