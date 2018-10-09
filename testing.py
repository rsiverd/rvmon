
##--------------------------------------------------------------------------##
## Examine how eclipse phases change with orbital ecc, W:
sys.stderr.write("-----------------------------------------------------\n")
npts_e, npts_w = 100, 100
e_list = np.linspace(0.0, 1.0, npts_e, endpoint=False)
w_list = np.linspace(0.0, 2.0*np.pi, npts_w, endpoint=False)
cphase = np.zeros((2, w_list.size, e_list.size), dtype='float32')


npoints = 10000
orb_phase = np.linspace(0.0, 1.0, npoints, endpoint=False)
mean_anom = 2.0 * np.pi * orb_phase
conj_true_anoms = [0.5*np.pi, 1.5*np.pi]    # anomalies at conjunction
#mean_anom = np.linspace(0.0, 2.0*np.pi, 100, endpoint=False)
for j,ecc in enumerate(e_list, 0):
    sys.stderr.write("Eccentricity %d of %d ...\n" % (j+1, e_list.size))
    true_anom = orbit.calc_true_anom(mean_anom, ecc)
    for i,w in enumerate(w_list, 0):
        LoS_anoms = true_anom + w 
        for k,anom in enumerate(conj_true_anoms):
            ecl_which = (qfo._fastRadianSep(LoS_anoms, anom)).argmin()
            #ecl_Manom = mean_anom[ecl_which]
            #cphase[k, j, i] = orb_phase[ecl_which]
            cphase[k, i, j] = mean_anom[ecl_which]


sys.stderr.write("-----------------------------------------------------\n")


XX, YY = np.meshgrid(e_list, w_list)
U0, V0 = np.cos(cphase[0]), np.sin(cphase[0])
U1, V1 = np.cos(cphase[1]), np.sin(cphase[1])

clf()
plt.quiver(XX, YY, U0, V0, pivot='mid')
xlabel('ecc')
ylabel('omega')

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Novel way to identify phase shift?

def vel_intersect(spec_RVs, model_phi, model_RVs, ftol=1e-6):
    """Look for phase shift using intersection of model and data."""
    results = []
    nextval = np.roll(model_RVs, -1)
    absdiff = np.abs(nextval - model_RVs)

    for trv in spec_RVs:   
        abs_sum = np.abs(nextval - trv) + np.abs(trv - model_RVs)
        btw_idx = np.where(np.isclose(abs_sum, absdiff, rtol=ftol))[0]
        isects = []
        for ix1 in btw_idx:
            ix2 = (ix1 + 1) % len(model_phi)
            phi_step = (model_phi[ix2] - model_phi[ix1]) % 1.0
            vel_step = (model_RVs[ix2] - model_RVs[ix1])
            stepfrac = (trv - model_RVs[ix1]) / vel_step
            newphase = (model_phi[ix1] + stepfrac * phi_step) % 1.0
            isects.append(round(newphase, 5))
            pass
        results.append(isects)
    return results

v_isects = vel_intersect(rad_vel, tphase, tcurve)
# %timeit blarg = shape_compare(tphase, tcurve, orbit_phase, rad_vel)
# 1000 loops, best of 3: 490 µs per loop

def test_orbit_shape(spec_phi, spec_RVs, model_phi, model_RVs, ftol=1e-6):
    """
    Compare orbit model to RV data points. Calculate phase offset
    and estimate quality of shape match.
    """
    v_isects = vel_intersect(spec_RVs, model_phi, model_RVs, ftol=ftol)

    offsets = (np.array(v_isects) - spec_phi[:, None]) % 1.0
    off_avg = np.average(offsets, axis=0)
    off_std = np.std(offsets, axis=0)
    bestidx = np.argmin(off_std)

    best_offset = off_avg[bestidx]
    best_stddev = off_std[bestidx]
    return (best_offset, best_stddev)

asdf = test_orbit_shape(orbit_phase, rad_vel, tphase, tcurve)

##--------------------------------------------------------------------------##
## Examine how eclipse phases change with orbital ecc, W:
sys.stderr.write("-----------------------------------------------------\n")
npts_e, npts_w = 20, 20
e_list = np.linspace(0.0, 1.0, npts_e, endpoint=False)
w_list = np.linspace(0.0, 2.0*np.pi, npts_w, endpoint=False)
#cphase = np.zeros((2, w_list.size, e_list.size), dtype='float32')


kmult = 1.001
npoints = 1000
v0 = qfo.elem_dict['v0']
KK  = kmult * qfo.elem_dict['K']
mod_phase = np.linspace(0.0, 1.0, npoints, endpoint=False)
mean_anom = 2.0 * np.pi * mod_phase
#conj_true_anoms = [0.5*np.pi, 1.5*np.pi]    # anomalies at conjunction
#mean_anom = np.linspace(0.0, 2.0*np.pi, 100, endpoint=False)


sys.stderr.write("Brute-force shape-fit test ... \n")
tik = time.time()
booya = qfo._fast_eccw_grid(nspec_phase, nspec_radvel, e_list, w_list, v0, KK)
#results = []
#for j,ecc in enumerate(e_list, 0):
#    sys.stderr.write("\rEccentricity %d of %d ... " % (j+1, e_list.size))
#    true_anom = orbit.calc_true_anom(mean_anom, ecc)
#    for i,w in enumerate(w_list, 0):
#        #LoS_anoms = true_anom + w 
#        mod_radvel = v0 + KK * orbit.calc_radvel(true_anom, ecc, w, norm=True)
#        fit_qual = qfo._test_orbit_shape(nspec_phase, nspec_radvel,
#                mod_phase, mod_radvel)
#        pshift, stddev = fit_qual
#        results.append((ecc, w, pshift, stddev))
#        pass
#    pass
sys.stderr.write("done.\n")

#results = np.array(results)
results = np.array(booya)
tok = time.time()
sys.stderr.write("Test completed in %.3f seconds.\n" % (tok-tik))

grid_eccen, grid_omega, grid_shift, grid_resid = results.T

best_idx = grid_resid.argmin()

#best_eccen = grid_eccen[best_idx]
#best_omega = grid_omega[best_idx]
#best_shift = grid_shift[best_idx]
#best_resid = grid_resid[best_idx]
best_eccen, best_omega, best_shift, best_resid = results[best_idx]

#rect_eccen = grid_eccen.reshape(npts_e, npts_w)
#rect_omega = grid_omega.reshape(npts_e, npts_w)

#rect_pshift = grid_pshift.reshape(npts_e, npts_w)
#rect_stddev = grid_stddev.reshape(npts_e, npts_w)

sys.stderr.write("Best parameters:\n")
sys.stderr.write("ecc, omega = %7.5f, %7.5f\n" % (best_eccen, best_omega))
sys.stderr.write("Phase shift: %7.5f\n" % best_shift)
sys.stderr.write("Best stddev: %7.5f\n" % best_resid)

##--------------------------------------------------------------------------##
## Circular average:
def circ_avg_phase(phase):
    y_avg = np.average(np.sin(2. * np.pi * phase))  # Cartesian Y avg
    x_avg = np.average(np.cos(2. * np.pi * phase))  # Cartesian X avg
    a_avg = np.arctan2(y_avg, x_avg)                # average angle (rad)
    p_avg = (a_avg / (2. * np.pi)) % 1.0            # decimal day average
    #tshift = 0.5 - t_avg                            # offset from midnight
    #sys.stderr.write("y_avg: %10.5f\n" % y_avg)
    #sys.stderr.write("x_avg: %10.5f\n" % x_avg)
    #sys.stderr.write("a_avg: %10.5f\n" % a_avg)
    #sys.stderr.write("p_avg: %10.5f\n" % p_avg)
    return p_avg

## Find best pair permutation:
def smartypants(offsets, mincontrast=0.2):
    avg_phase = circ_avg_phase(offsets.flatten())
    sys.stderr.write("avg_phase: %10.5f\n" % avg_phase)
    delta_phi = qfo._fastPhaseSep(offsets, good_phase)
    pcontrast = np.abs(np.diff(delta_phi, axis=1))
    good_vals = np.array([off[dphi.argmin()] for off,dphi,keep in \
            zip(offsets, delta_phi, (pcontrast >= 0.1)) if keep])
    cln_phase = circ_avg_phase(good_vals)
    sys.stderr.write("cln_phase: %10.5f\n" % cln_phase)
    return cln_phase




