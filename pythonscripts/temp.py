import numpy as np
import sys
import os
import healpy as hp
from numpy.linalg import pinv

tsteps = 259500
nbins = 300
Blm = np.full((tsteps, hp.sphtfunc.Alm.getsize(160)), 0.+0.j)
for file in os.listdir('blms_dipole/'):
	Blm[filenum*(tsteps//nbins):(filenum+1)*(tsteps//nbins), :] = np.loadtxt('blms_dipole/'+file, dtype=complex)
skycov = np.sum(Blm, axis=0)
np.savetxt('output/skycov.txt', skycov)

K = Blm.T
Cinv = K @ Blm
coupling = pinv(Cinv)

np.savetxt(f'output/coupling.txt', coupling)
