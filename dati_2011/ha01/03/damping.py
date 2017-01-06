from math import atan2, sqrt, pi
from random import random
fmt = r"%d & %d & %8.5f & %9.5f \\"

def modi():
    return 1+0.01*(random()-0.5)

freq = 8.43
omeg = 2 * pi * freq
ome2 = omeg * omeg
mass = 80000
stif = mass * ome2
pzer = 1200.0
dsta = pzer / stif
zeta = 0.038

print "m", mass/2
print "k", stif/2

for n, ffre in enumerate(( 40, 50, 60, 70)):
    beta = ffre / omeg
    dina = 1/sqrt((1-beta*beta)**2+(2*beta*zeta)**2)
    ampl = dsta*dina*1000000*modi()
    phas = atan2(2*beta*zeta,1-beta*beta)*modi()*180./pi
    print fmt % (n+1, ffre, ampl, phas)
