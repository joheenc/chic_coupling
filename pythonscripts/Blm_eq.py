import numpy as np
import healpy as hp
import os
Blm = []
numfiles = len(os.listdir('blms_eq'))
i = 1
for file in os.listdir('blms_eq'):
	Blm.append(np.loadtxt('blms_eq/'+file, dtype=complex))
	print(f'{i}/{numfiles}')
	i += 1
Blm = np.vstack(Blm)
skycov_approx = hp.sphtfunc.alm2map(np.sum(Blm, axis=0), nside=32)
np.savetxt('skycov_eq.txt', skycov_approx)
