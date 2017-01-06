from math import pi, sin, sqrt
"""Integrate the equation of motion as required in problem n.2 of 1st
homework 2010-11 using the linear acceleration algorithm."""

# how we format the data that will be computed step-by-step
fmt  = "%+10.4f " + 4*"%+10.4g " + "%+10.4f"

# problem data, given or computed in part 1 of the problem
m    = 35000.0   # mass of machine, kg
c    = 0.000     # assumed damping, Newton/(meter/second)
p_0  = 1000.0    # applied force, Newton
freq = 2*pi*5.0  # its circular freq., rad/second
t_0  = 6.0       # duration of transient, second
dur  = 12.0      # total duration, second
h    = 0.010     # time step, second
k    = 7971604.0 # design stiffness, Newton/meter

# define the load function, uses global values t_0, p_0, freq
def load(t):
    if t > t_0: return p_0 * sin(freq*t)
    else: return p_0*t/t_0 * sin(freq*t*t/t_0/2)

# linear acceleration coefficients
a_fac = c*h/2.0 + 3.0*m;   v_fac = 3.0*c + 6.0*m/h
stiff = k + (3.0*c + 6.0*m/h)/h

# initial state
t = 0
x, v, p = 0.0, 0.0, load(t)
a = (p - c*v - k*x) / m

# note that we know every thing at t=0, so that we can print it...
# print a header
print "#      t/s        x/m       fs/N",
print "    v/(m/s)  a/(m/s/s)       p/N"
# print the initial state
print fmt % (t, x, x*k, v, a, p)

# iterate with linear acceleration algorithm
while (t + h/2) < dur:
    # new t
    t  = t + h
    # find the increments
    dp = load(t) - p
    dx = (dp + a*a_fac + v*v_fac)/stiff
    dv = 3.0*dx/h - 3*v - a*h/2.0
    # update the state for the new t
    x, v, p = x + dx, v + dv, p + dp
    a  = (p - c*v - k*x) / m
    # print the current state
    print fmt % (t, x, x*k, v, a, p)
    
