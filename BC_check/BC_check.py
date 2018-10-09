
import sys
import astropy.time as astt
import numpy as np

test_file = 'datasets/spectra_KEBC09C05946.txt'
data = np.genfromtxt(test_file, dtype=None, names=True, delimiter=',')

t_start = astt.Time(data['obsdate'], scale='utc', format='isot')
halfexp = astt.TimeDelta(data['exptime'], format='sec') / 2.0
jd_mid  = (t_start + halfexp).jd

site_lat = -30.1673305556
site_lon = -70.8046611111   # East!
site_hei = 2201.0           # meters

star_ra_hrs = 15.92052182308321
star_ra_deg = 15.0 * star_ra_hrs
star_de_deg = 26.48410350592006

sys.stderr.write("Star coordinates: (%.7f , %.7f)\n"
        % (star_ra_deg, star_de_deg))

sys.stderr.write("BC calculator inputs:\n")
for stuff in zip(jd_mid, data['bcval']):
    sys.stderr.write("%.8f %12.6f\n" % stuff)

