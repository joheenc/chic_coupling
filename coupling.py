import numpy as np
import healpy as hp

def compute_orbit(orbdec0=0, orbdec1=0, axrotPerOrb=20, ntimes=1000, norbits=20):
    dslope = (orbdec1 - orbdec0) / ntimes
    phis = np.linspace(0., norbits*360., ntimes)
    return np.vstack([axrotPerOrb*phis, -phis, orbdec0+np.linspace(0., dslope, ntimes)]).T

def compute_B(beamfile, orbit, lmax=15):
    simbeam = hp.fitsfunc.read_map(beamfile)
    alms = hp.map2alm(simbeam, lmax=lmax)
    B = np.full((len(orbit), int((lmax+2)*(lmax+1)/2.0)), 0.+0.j)
    for rstep in range(len(orbit)):
        rot = hp.Rotator(rot=orbit[rstep], deg=True, eulertype='ZYX')
        B[rstep] = rot.rotate_alm(alms, inplace=False)
    return B
