import cubetools as ct
import sys

samp = ct.read_cube(open(sys.argv[1],'r'))
ct.freq_cutoff(samp,freq_cutoff=0.6)
ct.gaussian_averager(samp,sigma=6,nbr_dist=10,repeat=1)
ct.write_cube(samp,open("smoothed/"+sys.argv[1],"w"))
