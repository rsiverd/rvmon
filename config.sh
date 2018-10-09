## Vanderbilt paths:
shredder_hd1="/net/shredder.phy.vanderbilt.edu/hd1"
arch_root="${shredder_hd1}/siverd/arch_eng"
save_root="${shredder_hd1}/siverd/mayhem_proc"

## File and metadata lists:
recent_data="recent_NRES_spectra.txt"
recent_hdrs="recent_NRES_headers.txt"
recent_tars="recent_NRES_tarballs.txt"
recent_clean_data="recent_clean_NRES_spectra.txt"
recent_clean_hdrs="recent_clean_NRES_headers.txt"
recent_clean_tars="recent_clean_NRES_tarballs.txt"
all_objects="all_objects.txt"
object_list="object_list.txt"
blacklist="blacklist.txt"

## Where to put data points and plots:
plot_dir="orbit_plots"
meta_dir="datasets"
mkdir -p $meta_dir $plot_dir

