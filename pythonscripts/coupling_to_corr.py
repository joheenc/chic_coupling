import numpy as np
import healpy as hp
from coupling import compute_corr
import sys

lmax = int(sys.argv[1])
coupling = np.loadtxt('output/coupling.txt', dtype=complex)
corrs = compute_corr(coupling, lmax=lmax)
np.savetxt('output/corrs_a10.txt', corrs)
