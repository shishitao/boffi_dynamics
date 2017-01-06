from math import sqrt, sin, cos, exp, pi

# data
k = 50000.0 ; m = 1200.00 ; z = 0.05
# p = 1900 + 200 t/t0
def pt(t,t0=0.5): return 1900.0+400.0*t if t<=t0 else 0.0

# time step
h = 0.02                       # time step

# derived system parameters
w2 = k/m ; wn = sqrt(w2) ; wd = wn * sqrt(1-z*z) ; c  = 2*m*wn*z

# magical constants for Linear Acceleration algorithm
k_hat = k + 3*c/h + 6*m/h/h
c_hat = 3*c + 6*m/h
m_hat = c*h/2 + 3*m
dx2dv = 3.0/h ; v2dv = 3.0 ; a2dv = h/2.0

# initial conditions
x = 0 ; v = 0 ; t = 0 ; p = pt(t)
a = p - c*v - k*x ; a = a/m

# print values, advance solution
while t<2+h/2:
    
    print (4*" %+014.8f") % (t, x*1000, v*1000, p)

    t1 = t+h
    p1 = pt(t1)
    dp = p1-p
    dp_hat = dp + c_hat*v + m_hat*a

    dx = dp_hat/k_hat
    dv = dx*dx2dv - v*v2dv - a*a2dv

    x = x + dx ; v = v +dv ; p = p1 ; t = t1
    a = p - c*v - k*x ; a = a/m
