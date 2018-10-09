## Vanderbilt paths:
shredder_hd1="/net/shredder.phy.vanderbilt.edu/hd1"
arch_root="${shredder_hd1}/siverd/arch_eng"
save_root="${shredder_hd1}/siverd/mayhem_proc"

## NRES sites and obs-months to use:
nres_sites=( lsc elp cpt tlv )
#month_list=( 2018{01..08} )
month_list=( '??????' )

## File and metadata lists:
info_dir="meta_info"
mkdir -p $info_dir
recent_all_data="${info_dir}/recent_all_NRES_spectra.txt"
recent_all_img_hdrs="${info_dir}/recent_all_NRES_headers.txt"
recent_all_tarballs="${info_dir}/recent_all_NRES_tarballs.txt"
recent_cln_data="${info_dir}/recent_cln_NRES_spectra.txt"
recent_cln_img_hdrs="${info_dir}/recent_cln_NRES_headers.txt"
recent_cln_tarballs="${info_dir}/recent_cln_NRES_tarballs.txt"
all_objects="all_objects.txt"
object_list="object_list.txt"
blacklist="blacklist.txt"

## Where to put data points and plots:
plot_dir="orbit_plots"
meta_dir="datasets"
mkdir -p $meta_dir $plot_dir

