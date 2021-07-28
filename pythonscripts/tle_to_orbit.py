from coupling import orbit_from_tle
import numpy as np
import sys

tlefile = sys.argv[1]
tle = []
with open(tlefile, 'r') as infile:
	for line in infile.readlines():
		tle.append(line.strip())
tstart = float(sys.argv[2]) #in MJD
tstop = float(sys.argv[3])
bin = int(sys.argv[4])
nbins = int(sys.argv[5])
tstarts = np.linspace(tstart, tstop, nbins+1)
ras, decs = orbit_from_tle(tle[0], tle[1], tstarts[bin], tstarts[bin+1])
np.savetxt(f'orbitfiles/orbit_bin_{bin}.txt', np.vstack([ras, decs]))
