#!/usr/bin/python
from matplotlib.pyplot import plot,title,annotate,gcf
from numpy import array
import sys

def readgeo(inpf):
  dat = []

  cur = inpf.readline().split()
  while cur != 'ATOM X/A Y/B Z/C'.split():
    cur = inpf.readline().split()
  inpf.readline()
  cur = inpf.readline().split()

  while cur != []:
    cur.pop(3)
    cur.pop(1)
    dat.append(cur)
    cur = inpf.readline().split()

  return array(dat).astype(float).T


atis, atyps, ax, ay, az = readgeo(open(sys.argv[1],'r'))
ax%=1; ay%=1; az%=1

nfe = list(atyps).count(226)
points = list(zip(ax[0:nfe],ay[0:nfe]))
slabs  = atis[0:nfe].astype(int)

plot(ax[0:nfe],ay[0:nfe],'k.')
title(sys.argv[1])

for i in range(len(slabs)):
  annotate(slabs[i], points[i],size='small')

gcf().savefig("geo_display.pdf")

