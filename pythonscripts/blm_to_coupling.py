import numpy as np
import sys
import os
import healpy as hp
from numpy.linalg import pinv
from pygdsm import GlobalSkyModel
from coupling import gal_to_dipole

tsteps = 259500
nbins = 300
nside = 64
lmax = int(sys.argv[1])

Blm = np.full((tsteps, hp.sphtfunc.Alm.getsize(lmax)), 0.+0.j)
filenum = 0
for file in os.listdir('blms_dipole/'):
	Blm[filenum*(tsteps//nbins):(filenum+1)*(tsteps//nbins), :] = np.loadtxt('blms_dipole/'+file, dtype=complex)
	filenum += 1
skycov = np.sum(Blm, axis=0)
np.savetxt('output/skycov.txt', skycov)

alm = hp.map2alm(gal_to_dipole(hp.ud_grade(GlobalSkyModel(freq_unit='MHz').generate(150), nside)), lmax=lmax)
np.savetxt('output/sky.txt', alm)
P = Blm @ alm #unmapped output
np.savetxt('output/P.txt', P)

K = Blm.T
coupling = K @ Blm
np.savetxt('output/coupling.txt', coupling)

Imap = coupling @ alm #mapped output
np.savetxt('output/Imap.txt', Imap)
