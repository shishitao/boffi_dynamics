from math import *

def resp_elas(m,c,k, cC,cS,w, F, x0,v0):
    # response (disp+vel) of a linear SDOF subjected to a harmonic load
    # use:
    # x,v=resp_elas(...) ; vel_35=v(3.5)
    # NB: cC and cS are FORCES, so that p(t) = cC cos(w t) + cS sin(w t)
    #     cC and cS may also be equal to zero, so that we have free response
    wn2=k/m ; wn=sqrt(wn2) ; beta=w/wn
    z=c/(2*m*wn)
    wd=wn*sqrt(1-z*z)
    # csi(t) = R sin(w t) + S cos(w t) + D
    det=(1.-beta**2)**2+(2*beta*z)**2
    R=((1-beta**2)*cS + (2*beta*z)*cC)/det/k
    S=((1-beta**2)*cC - (2*beta*z)*cS)/det/k
    D=F/k
    # x(t) = exp(-z wn t)(A cos(wd t) + B sin(wd t)) + csi(t) + D
    # v(t) = exp(-z wn t) *
    #        [ + wd * (B cos(wd t) - A sin(wd t)) -z wn * (A cos(wd t)+B sin(wd t)) ] +
    #            w R cos(w t) - w S sin(w t)
    # x(0) = 1 * ( A*1 + B*0 ) + R*0 + S*1 + D = x0
    A=x0-S-D
    # v(0) = wd B - z wn A + w R = v0
    B=(v0+z*wn*A-w*R)/wd
    def x(t):
        return exp(-z*wn*t)*(A*cos(wd*t)+B*sin(wd*t))+R*sin(w*t)+S*cos(w*t)+D
    def v(t):
        return (-z*wn*exp(-z*wn*t)*(A*cos(wd*t)+B*sin(wd*t))
                  +wd*exp(-z*wn*t)*(B*cos(wd*t)-A*sin(wd*t))
                  +w*(R*cos(w*t)-S*sin(w*t)))
    return x,v

def resp_yeld(m,c, cC,cS,w, F, x0,v0):
    # csi(t) = R sin(w t) + S cos(w t) + \alpha t
    alpha=F/c
    det=w**2*(c**2+w**2*m**2)
    R=(+w*c*cC-w*w*m*cS)/det
    S=(-w*c*cS-w*w*m*cC)/det
    # x(t) = A exp(-c t/m) + B + R sin(w t) + S cos(w t) + alpha t
    # v(t) = - c A/m exp(-c t/m) + w R cos(w t) - w S sin(w t) + alpha
    #
    # v(0) = -c A / m + w R + alpha = v0
    A=m*(w*R+alpha-v0)/c
    # x(0) = A + B + S = x0
    B=x0-A-S
    def x(t):
        return A*exp(-c*t/m)+B+R*sin(w*t)+S*cos(w*t)+alpha*t
    def v(t):
        return -c*A*exp(-c*t/m)/m+w*R*cos(w*t)-w*S*sin(w*t)+alpha
    return x,v

def shift(f,t0):
    def f1(t):
        return f(t-t0)
    return f1

def bisect(f,t,x0,x1,count):
    h=(x0+x1)/2.0
    fh=f(h)-t
    if abs(fh)<1e-8 : return h
    f0=f(x0)-t
    f1=f(x1)
    if f0*fh>0 : return bisect(f,t,h,x1,count+1)
    return bisect(f,t,x0,h,count+1)

def printrange(f,origin,duration,nsteps):
    dtau=duration/nsteps
    for i in range(nsteps+1):
        tau=dtau*i
        t=origin+tau
        print t, f(tau)

# SDOF characteristics
mass=1000. # kg
k=40000.   # N/m
zeta=0.03  # damping ratio
damp=2*zeta*mass*sqrt(k/mass)
fy=2500.   # N

xy=fy/k    # m
 
# load characteristics, p(t) = Po sin(pi t/t1) + 0 cos(pi t/t1), 0 <= t <= t1
#                            = 0                                  otherwise.
t1=0.3     # s
w=pi/t1    # rad/s
Po=6000.   # N

# initial conditions
x0=0.0     # m
v0=0.0     # m/s
# initially elastic, compute the response functions
x_next,v_next=resp_elas(mass,damp,k, 0.0,Po,w, 0.0, x0,v0)

# find the time when the spring yelds --- it happens that's LESS than t1
t_yeld=bisect(x_next,xy,0.0,t1,1)

# scan the range 0 ___ t_yeld and print
dtau=(t_yeld-0.0)/100.
for i in range (0,101):
    tau=dtau*i
    t=tau+0.0
    print t, x_next(tau)

# we're yelded, find the new initial conditions
x0=x_next(tau)
v0=v_next(tau)
# the coefficients of the loading function with origin in t_yeld
cS=cos(w*tau)*Po
cC=sin(w*tau)*Po
# compute the response functionf in the yelding phase, with sin load
x_next,v_next=resp_yeld(mass,damp,    cC,cS,w, -fy, x0,v0)

dt=(t1-t_yeld)/100.
for i in range(101):
    tau=i*dt
    t=tau+t_yeld
    print t, x_next(tau)

# after t1, p(t)=0.0
cS=0.0 ; cC=0.0

# initial conditions for yelded phase, but p=0.0
x0=x_next(tau)
v0=v_next(tau)
# compute the response functions in yeld, with no load
x_next,v_next=resp_yeld(mass,damp,    cC,cS,w, -fy, x0,v0)

# find the time for which the velocity is zero
t2=t1+bisect( v_next, 0.0, 0, 0.3, 1.0)

dt=(t2-t1)/100
for i in range(101):
    tau=i*dt
    t=tau+t1
    print t, x_next(tau)

# now the velocity is 0.0, going back to elastic behaviour
x0=x_next(tau) ; v0=0.0
# the constant force models the permanent displacements
x_next,v_next=resp_elas(mass,damp,k, 0.0,0.0,w, k*x0-fy, x0,v0)

dt=(4.0-t2)/200
for i in range(201):
    tau=i*dt
    t=tau+t2
    print t, x_next(tau)
