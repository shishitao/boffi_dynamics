from scipy import *

# solves part 1 of problem #3, 1st home work 2010-11

a = matrix("1 -1600;1 -2500;1 -3600;1 -4900")
cosines = cos(array((7.58258,33.33505,163.21210,171.69968))*pi/180)
rho     = array((12.39062,41.09556,18.07490,7.11246))/1E6
p_0     = 600

ata = a.T*a
print ata
b   = matrix(p_0*cosines/rho).T
print b/1E6
print a.T*b/1E6
x   = ata.I*a.T*b
print x

# for a single variable (case of \zeta estimation) the least square
# approximation is the arithmetic mean of the estimates

print "3.8%"
