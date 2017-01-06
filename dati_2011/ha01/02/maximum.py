from math import pi, sin, sqrt
"""Integrate the eom for problem 2 of 1st homework 2010-11, using
different values  of the external force transient duration and of the
damping ratio, printing only the peak value of the transmitted force."""

# Data
m    = 35000.0   # kg

p_0  = 1000.0    # Newton
p_m  =  300.0    # Newton

freq = 5.0       # Hertz
cfre = freq*2*pi # rad/second
dur  = 12.0      # second
h    = 0.010     # second

# cycle over damping ratio, it is  0.000 <= z <= 0.200
for i in range(101):
    z = i*0.002
    # compute the maximum stiffness allowed, see solutions'paper
    # for the derivation
    k = 1E6*(sqrt((182*z**2 + 9)**2 + 819) - (182*z**2 + 9) )*pi**2 / 26
    # we know dampin ratio, we need damping coefficient
    w = sqrt(k/m) # rad/second
    c = 2*z*w*m   # Newton/(meter/second)

    # linear acceleration coefficients for this value of z
    a_fac = c*h/2.0 + 3.0*m
    v_fac = 3.0*c + 6.0*m/h
    stiff = k + (3.0*c + 6.0*m/h)/h
        
    # cycle over external force transient duration t_0
    # it is 2.00s <= t_0 <= 8.00s
    for j in range(101):
        t_0 = 2+j*0.06

        def load(t):
            if t < t_0:
                return p_0*t/t_0*sin(cfre*t*t/t_0/2.0)
            else:
                return p_0 * sin(cfre*t)

        # initial state
        x, v, p = 0.0, 0.0, load(0.0)
        a       = (p - c*v - k*x) / m
        t       = 0.0
        f_max   = 0.0
        
        while True:
            t  = t+h
            if t-h/2. > dur: break
            dp = load(t)-p
            dx = (dp+a*a_fac+v*v_fac)/stiff
            dv = 3.0*dx/h - 3*v - a*h/2.0
            x  = x+dx
            v  = v+dv
            p  = p+dp
            a  = (p - c*v - k*x) / m
            f_max = max(f_max, abs(x*k+v*c))

        print "%5.3f %6.4f %8.2f" % (z, t_0, f_max)

    # the plotting program to be used wants a blank line between
    # blocks of data with the same abscissa, make it content
    print 

# the two loops, inner on t_0 outer on z, are done
