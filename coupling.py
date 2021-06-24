import numpy as np
import healpy as hp

def compute_orbit(orbdec0=0, orbdec1=0, axrotPerOrb=20, ntimes=1000, norbits=20):
	dslope = (orbdec1 - orbdec0) / ntimes
	phis = np.linspace(0., norbits*360., ntimes)
	return np.vstack([axrotPerOrb*phis, -phis, orbdec0+np.linspace(0., dslope, ntimes)]).T

def compute_B(beamfile, orbdec0=0, orbdec1=0, axrotPerOrb=20, ntimes=1000, norbits=20, lmax=15): #compute the beam matrix for a particular frequency
	simbeam = hp.fitsfunc.read_map(beamfile)
	alms = hp.map2alm(simbeam, lmax=lmax)
	dslope = (orbdec1 - orbdec0) / ntimes
	phis = np.linspace(0., norbits*360., ntimes)
	B = np.zeros((ntimes, int((lmax+2)*(lmax+1)/2.0))) #B matrix stores alms (columns) at each t step (rows)
	for rstep in range(ntimes):
		rot = hp.Rotator(rot = [axrotPerOrb*phis[rstep], -phis[rstep], orbdec0+rstep*dslope], deg=True, eulertype='ZYX')
		B[rstep] = rot.rotate_alm(alms, inplace=False))
	return B
