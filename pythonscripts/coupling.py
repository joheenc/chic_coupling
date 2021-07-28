import numpy as np
import healpy as hp
from scipy.ndimage import rotate
from sgp4.api import Satrec, SGP4_ERRORS
from astropy.coordinates import SkyCoord, ITRS, TEME, CartesianDifferential, CartesianRepresentation
from astropy.time import Time
from astropy import units as u

def point_beam(theta, phi, beamthetastd=0.3, beamphistd=0.3, nres=400, nside=32, invert=True, savefile=None):
    inv = -1 if invert else 1
    thetagrid = np.hstack(np.repeat(np.linspace(0, np.pi, nres), nres).reshape(nres, nres))
    phigrid = np.hstack(np.repeat(np.linspace(-np.pi, np.pi, nres), nres).reshape(nres, nres).T)
    thetadist = thetagrid-theta
    phidist = np.minimum(np.minimum(np.abs(phigrid-inv*phi), np.abs(phigrid+2*np.pi-inv*phi)), \
                         np.abs(phigrid-2*np.pi-inv*phi)) #periodic boundary condition
    beam = 1/(2*np.pi*beamthetastd*beamphistd) * np.exp(-0.5*(thetadist**2/beamthetastd**2 + \
                                                              phidist**2/beamphistd**2))
    indices = hp.ang2pix(nside, thetagrid, phigrid, nest=False)
    hpxmap = np.zeros(hp.nside2npix(nside), dtype=np.float)
    hpxmap[indices] += beam
    if savefile != None:
        hp.fitsfunc.write_map(savefile, hpxmap)
    return beam, indices, hpxmap/np.max(hpxmap)

def rotate_beam(beam, indices, hpxmap, theta1, phi1, theta2, phi2, psi=0, nres=400, nside=32, onlybeam=False):
    thetas = np.linspace(0, np.pi, nres)
    phis = np.hstack(np.repeat(np.linspace(-np.pi, np.pi, nres), nres).reshape(nres, nres).T)
    theta1idx = np.argmin(np.abs(thetas-theta1))
    phi1idx = np.argmin(np.abs(phis-phi1))
    if psi != 0:
        centeredbeam = rotate_beam(beam, indices, hpxmap, theta1, phi1, np.pi/2, 0, psi=0, nres=nres, onlybeam=True)
        beamgrid = centeredbeam.reshape(nres, nres)
        beamgrid = rotate(beamgrid, angle=psi*180./np.pi, reshape=False) #from scipy.ndimage
        centeredbeam = beamgrid.flatten()
        beam = rotate_beam(centeredbeam, indices, hpxmap, np.pi/2, 0, theta1, phi1, psi=0, nres=nres, onlybeam=True)
    theta2idx = np.argmin(np.abs(thetas-theta2))
    phi2idx = np.argmin(np.abs(phis-phi2))
    beam = np.roll(beam, (theta2idx-theta1idx)*nres)
    beam = np.roll(beam, (-phi2idx+phi1idx))
    
    if onlybeam:
        return beam
    hpxmap.fill(0)
    hpxmap[indices] += beam
    return beam, hpxmap/np.max(hpxmap)

#returns an orbit path in (RA, dec) coordinates given Two-Line Element orbit data and start/stop times in MJD
def orbit_from_tle(tle1, tle2, starttime, stoptime, tstep=6.944e-4):
    satellite = Satrec.twoline2rv(tle1, tle2)
    mjds = np.arange(starttime, stoptime, tstep)
    ras = np.zeros(len(mjds))
    decs = np.zeros(len(mjds))

    step = 0
    for mjd in mjds:
        t = Time(mjd, format='mjd')
        error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2)  # in km and km/s
        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])
            
        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
        itrs = teme.transform_to(ITRS(obstime=t))
        decs[step] = itrs.earth_location.lat.value
#         lons[step] = itrs.earth_location.lon.value % 360
        ras[step] = (itrs.earth_location.lon.value % 360 + 360*(mjd%1)) % 360
        step += 1
    return ras, decs
    

def sky_coverage(lons, lats, psis, coordsys='galactic', beamthetastd=0.3, beamphistd=0.3, nres=400, nside=32, progress=0, display=False):
    skycov = np.zeros((lons.shape[0], hp.nside2npix(nside)))
    if coordsys in ['celestial', 'equatorial', 'c', 'C']: #conversion to galactic coords
        for i in range(lons.shape[0]):
            c = SkyCoord(lons[i], lats[i], frame='icrs', unit='deg')
            lons[i] = c.galactic.l.degree
            lats[i] = c.galactic.b.degree

    phis = (lons%360-180)*np.pi/180
    thetas = (90-lats)*np.pi/180
    beam, indices, beammap = point_beam(thetas[0], phis[0],\
                                        beamthetastd=beamthetastd, beamphistd=beamphistd, nres=nres, nside=nside)
    skycov[0] = beammap
    for i in range(1, lons.shape[0], 1):
        beam, beammap = rotate_beam(beam, indices, beammap, thetas[i-1], phis[i-1],\
                              thetas[i], phis[i], psi=psis[i], nres=nres, nside=nside)
        if display:
            hp.mollview(beammap)
        skycov[i] = beammap
        if progress != 0 and i%progress==0:
            print(i)
    return skycov

def gal_to_dipole(hpxmap): #convert galactic coordinates to "dipole" coordinates
    phi_cmb_dipole = 263.85 
    theta_cmb_dipole = (90-48.25)
    rot_dipole = hp.Rotator([phi_cmb_dipole, theta_cmb_dipole, 0], eulertype='ZXZ')
    return rot_dipole.rotate_map_pixel(hpxmap)

def compute_B(skycov, lmax=15):
    blm = np.array([hp.map2alm(rot_dipole.rotate_map_pixel(step), lmax=lmax) for step in skycov])
    return blm

def compute_corr(B, lref=1, mref=0, lmax=15):
    a10 = B[:, hp.Alm.getidx(lmax, lref, mref)]
    corrs = np.full(hp.Alm.getsize(lmax), 0.+0.j)
    for l in range(lmax):
        for m in range(l+1):
            alm = B[:, hp.Alm.getidx(lmax, l, m)]
            corrs[hp.Alm.getidx(lmax, l, m)] = np.correlate(a10, alm)
    return corrs
