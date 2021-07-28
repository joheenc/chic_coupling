from coupling import *
import sys

bin = int(sys.argv[1])
lmax = int(sys.argv[2])
radecs = np.loadtxt(f'orbitfiles/orbit_bin_{bin}.txt')
psis = np.linspace(0, 0, radecs.shape[1]) % 2*np.pi
skycov = sky_coverage(radecs[0, :], radecs[1, :], psis, beamthetastd=0.2, beamphistd=0.4, coordsys='c', nside=64)
np.savetxt(f'skycovfiles/skycov_bin_{bin}.txt', skycov)
#skycov_dipole = gal_to_dipole(np.sum(skycov, axis=0))

#blm_eq = np.array([hp.map2alm(bstep, lmax=lmax) for bstep in skycov])
#np.savetxt(f'blms_eq/blm_eq_{bin}.txt', blm_eq)

blm_dipole = np.array([hp.map2alm(gal_to_dipole(bstep), lmax=lmax) for bstep in skycov])
np.savetxt(f'blms_dipole/blm_dipole_{bin}.txt', blm_dipole)
