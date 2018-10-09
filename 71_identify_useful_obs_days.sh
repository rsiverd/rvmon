#!/bin/bash
#
# Identify DAY-OBS relevant to targets of interest. These should get priority
# when rebuilding/reprocessing data.
#
# Rob Siverd
# Created:      2018-08-07
# Last updated: 2018-08-07
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.10"
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
baz="$tmp_dir/baz_$$.txt"
qux="$tmp_dir/qux_$$.txt"
jnk="$foo $bar $baz $qux"  # working copy
def_jnk="$jnk"             # original set
dir_cleanup='(echo -e "\nAutomatic clean up ... " ; cmde "rm -vrf $tmp_dir")'
jnk_cleanup='for X in $jnk ; do [ -f $X ] && cmde "rm -vf $X" ; done'
trap "$jnk_cleanup >&2" EXIT
##trap '[ -d $tmp_dir ] && cmde "rm -vrf $tmp_dir"' EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup >&2" EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup >&2; $jnk_cleanup >&2" EXIT
#trap 'oops=$? ; echo ; exit $oops' HUP INT TERM

## Required programs:
declare -a need_exec
need_exec+=( awk cat FuncDef sed tr )
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
#target_list="$1"

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

nres_sites=( "lsc" "elp" "cpt" "tlv" )
nres_cameras=( nres{01..04} )
header_data="recent_NRES_headers.txt"
target_list="useful_objects.txt"
[ -f $target_list ] || PauseAbort "Can't find file: $target_list"
[ -f $header_data ] || PauseAbort "Can't find file: $header_data"

## Load targets from file:
targets=( `cat $target_list` )
total=${#targets[*]}

## File modification age check:
seconds_since_modification () {
   last_change_sec=`stat -c'%Z' $1`
   current_unix_sec=$(date -u +%s)
   mod_age_seconds=$(( $current_unix_sec - $last_change_sec ))
   echo $mod_age_seconds
}

## Check age of header_data file:
max_age_sec=86400
updated_ago_sec=$(seconds_since_modification $header_data)
echo "updated_ago_sec: $updated_ago_sec"
if [ $updated_ago_sec -ge $max_age_sec ]; then
   recho "Header metadata stale and need refresh ...\n"
   cmde "bash ./EVERYDAY.sh"
else
   gecho "Header data still current!\n"
fi

## Look for matches:
rm $foo 2>/dev/null
for item in ${targets[*]}; do
   grep $item $header_data >> $foo
done
cmde "wc -l $foo"
echo

## Create a folder for lists of useful obs days:
dayobs_dir="rjs_dayobs"
cmde "mkdir -p $dayobs_dir" || exit $?
echo

## Get site,day-obs pairs that warrant reprocessing:
#rm $bar 2>/dev/null
#proc_lists=()
for nrcam in ${nres_cameras[*]}; do
   days_list="$dayobs_dir/${nrcam}.txt"
   yecho "Populating $days_list ... "
   #echo "nrcam: $nrcam"
   grep $nrcam $foo | awk '{ print $1 }' | sed 's|^.*/||g' > $baz
   nhits=$(cut -d- -f3 $baz | sort -u | tee $qux | wc -l)
   echo "nhits: $nhits"
   if [ $nhits -gt 0 ]; then
      #echo "$nrcam $(cat $qux | tr '\n' ' ')" >> $bar
      cmde "mv -f $qux $days_list" || exit $?
      #proc_lists+=( $days_list )
   fi
done
#cmde "mv -f $bar $proc_list" || exit $?
#cmde "cat $bar"
#exit 0

## Save a copy of reprocessing recommendations to mayhem folder:
yecho "\nCopying day-obs lists to pipeline folder ...\n"
mayhem_dayobs_dir="$HOME/NRES/mayhem/rjs_obsdays"
cmde "rsync -av $dayobs_dir/ $mayhem_dayobs_dir/"
#yecho "Saving DAY-OBS list copies to mayhem folder ...\n"
#cmde "mkdir -p $mayhem_dayobs_dir"                 || exit $?
#for item in ${proc_lists[*]}; do
#   save_days="$mayhem_dayobs_dir/$(basename $item)"
#   cmde "cp -f $item $qux"                         || exit $?
#   cmde "mv -f $qux $save_days"                    || exit $?
#done


##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
[ -f $baz ] && rm -f $baz
[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (71_identify_useful_obs_days.sh):
#---------------------------------------------------------------------
#
#  2018-08-07:
#     -- Increased script_version to 0.10.
#     -- First created 71_identify_useful_obs_days.sh.
#
