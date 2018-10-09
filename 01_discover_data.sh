#!/bin/bash
#
# Identify raw (e00) and processed (e91) NRES products to be used by
# subsequent monitoring scripts. File location information shall be
# picked up from config.sh.
#
# Rob Siverd
# Created:      2018-04-17
# Last updated: 2018-08-08
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.05"
this_prog="${0##*/}"
#shopt -s nullglob
# Propagate errors through pipelines: set -o pipefail
# Exit if uninitialized variable used (set -u): set -o nounset
# Exit in case of nonzero status (set -e): set -o errexit

## Program options:
#save_file=""
#shuffle=0
#confirmed=0

## Standard scratch files/dirs:
tmp_name="$(date +%Y%m%d.%H%M%S).$$.$(whoami)"
tmp_root="/tmp"
[ -d /dev/shm ] && [ -w /dev/shm ] && tmp_root="/dev/shm"
tmp_dir="$tmp_root"
#tmp_dir="$tmp_root/$tmp_name"
#mkdir -p $tmp_dir
foo="$tmp_dir/foo_$$.txt"
bar="$tmp_dir/bar_$$.txt"
baz="$tmp_dir/baz_$$.fits"
qux="$tmp_dir/qux_$$.fits"
jnk="$foo $bar $baz $qux"  # working copy
def_jnk="$jnk"             # original set
dir_cleanup='(echo -e "\nAutomatic clean up ... " ; cmde "rm -vrf $tmp_dir")'
jnk_cleanup='for X in $jnk ; do [ -f $X ] && cmde "rm -vf $X" ; done'
trap "$jnk_cleanup" EXIT
##trap '[ -d $tmp_dir ] && cmde "rm -vrf $tmp_dir"' EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup" EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup ; $jnk_cleanup" EXIT
#trap 'oops=$? ; echo ; exit $oops' HUP INT TERM

## Required programs:
declare -a need_exec
need_exec+=( awk cat FuncDef iltk sed tr )
#need_exec+=( shuf shuffle sort ) # for randomization
for need in ${need_exec[*]}; do
   if ! ( /usr/bin/which $need >& /dev/null ); then
      echo "Error:  can't find '$need' in PATH !!" >&2
      exit 1
   fi
done

## Helper function definitions:
fd_args="--argchk --colors --cmde --echo"
#fd_args+=" --Critical"
fd_args+=" --rowwrite"
#fd_args+=" --timers"
fd_args+=" --warnings"
FuncDef $fd_args >/dev/null || exit $?
eval "$(FuncDef $fd_args)"  || exit $?

## Check for arguments:
usage () { 
   Recho "\nSyntax: $this_prog --START\n\n"
   #Recho "\nSyntax: $this_prog arg1\n\n"
}
if [ "$1" != "--START" ]; then
#if [ -z "$1" ]; then
   usage >&2
   exit 1
fi

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

## Load config:
cmde "source config.sh" || exit $?
touch $blacklist

## Make sure arch_root is defined and present:
[ -z "$arch_root" ] && ErrorAbort "Data storage path 'arch_root' not defined!"
[ -d $arch_root ] || ErrorAbort "Can't find directory: $arch_root"

## Ensure nres_sites is defined and has non-zero length:
[ ${#nres_sites[*]} -eq 0 ] && ErrorAbort "Empty nres_sites in config.sh!"
[ ${#month_list[*]} -eq 0 ] && ErrorAbort "Empty month_list in config.sh!"

##--------------------------------------------------------------------------##

## Find files, extract metadata:
keys="DATE-OBS EXPTIME BLKUID MOLUID OBJECTS"
arch_dirs=()
yecho "Searching "
rm $foo $bar 2>/dev/null
for site in ${nres_sites[*]}; do
   yecho "$site ... "
   arch_site="$arch_root/$site"
   for month in ${month_list[*]}; do
      find $arch_site/nres0?/${month}?? -type f -name "*e00.fits.fz" 2>/dev/null >> $foo
      find $arch_site/nres0?/${month}?? -type f -name "*e91.tar.gz"  2>/dev/null >> $bar
   done
done
gecho "done.\n"

## Store completed lists:
cmde "mv -f $foo $recent_all_data" || exit $?
cmde "mv -f $bar $recent_all_tarballs" || exit $?
cmde "wc -l $recent_all_data $recent_all_tarballs"
#exit
#find ${arch_dirs[*]} -type f -name "*e00.fits.fz" > $recent_data
#find ${arch_dirs[*]} -type f -name "*e91.tar.gz"  > $recent_tarballs

## Collect header keywords from image files:
cmde "imhget $keys --progress -l $recent_all_data -o $foo"      || exit $?
cmde "mv -f $foo $recent_all_hdrs"                              || exit $?

## Purge black-listed files:
cmde "iltk -L $recent_all_data -e $blacklist -o $recent_cln_data" || exit $?
cmde "iltk -L $recent_all_tarballs -e $blacklist -o $recent_cln_tarballs" || exit $?
cmde "iltk -L $recent_all_img_hdrs -e $blacklist -o $recent_cln_img_hdrs" || exit $?

## Extract objects:
cut -d' ' -f6- ${recent_cln_hdrs} | sed 's/\&thar.*$//' \
   | sort -u > $all_objects

##--------------------------------------------------------------------------##

## select KELT objects:
grep KEB ${recent_cln_hdrs} | awk '{ print $6 }' | cut -d\& -f1 \
   | sort -u > $object_list
cat <(grep KEB $all_objects) <(grep KELT $all_objects) \
    <(grep Jupiter $all_objects) <(grep rocyon $all_objects) \
    <(grep WASP $all_objects) <(grep HAT $all_objects) \
    <(grep '^KS' $all_objects) | awk '{ print $1 }' | sort -u > $object_list


##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
#[ -f $baz ] && rm -f $baz
#[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (01_discover_data.sh):
#---------------------------------------------------------------------
#
#  2018-08-08:
#     -- Increased script_version to 0.04.
#     -- Now mine headers into temporary files to avoid leaving headers
#           list in bad state.
#     -- Search for useful spectra extended to current date and now
#           includes nres03 and nres04 (CPT and TLV).
#
#  2018-06-11:
#     -- Increased script_version to 0.03.
#     -- Added WASP and HAT targets to the list of things processed.
#
#  2018-05-07:
#     -- Increased script_version to 0.02.
#     -- Now pull from 'clean' data files (blacklisted files removed).
#
#  2018-04-17:
#     -- Increased script_version to 0.01.
#     -- First created 01_discover_data.sh.
#
