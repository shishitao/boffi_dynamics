"""Solves exercise #5, 1st homework""" 
from scipy import mat

# everything is adimensional

# the dimensional coefficients attached to each expression are
# reminded in the comments

M=mat('4 1 0;1 4 1;0 1 2')/6.  # m
K=mat('1 0 0;0 2 0;0 0 3')*1.  # k
v=mat('1;1;1')                 # Z

f=M*v                # w^2 m Z
u=K.I*f              # w^2 m Z / k

V_0 = (v.T*K*v)[0,0] # Z^2 k
T_0 = (v.T*M*v)[0,0] # w^2 Z^2 m
V_1 = (f.T*u)[0,0]   # w^2 m Z w^2 m Z /k = w^4 m^2 Z / k
T_1 = (u.T*M*u)[0,0] # (w^3 m Z / k) m (w^3 m Z / k) = w^6 m^3 Z^2 / k^2

print "R_oo =", V_0/T_0, "k/m"
print "R_oi =", T_0/V_1, "k/m"
print "R_ii =", V_1/T_1, "k/m"
