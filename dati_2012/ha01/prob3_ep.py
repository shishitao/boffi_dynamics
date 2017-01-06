from math import pi, sqrt, sin
fmt="%6.3f %9.4f    %+8.4f %+8.4f %+8.4f     %+8.2e %+8.4f %+8.2e"

def make_p(p0,p1,t0,t1):
  """make_p(p0,p1,t0,t1) returns a zero, linear, zero 1/2 load function, p(t)"""
  dp = p1-p0 ; dt = t1-t0 ; a = dp/dt ; b = p0 - a*t0
  def p(t):
    if t0<=t<=t1:
      return a*t+b
    else:
      return 0.0
  return p

def make_kt(k,fy):
  "make_kt(k,fy) returns a function kt(u,v,up) returning kt, up"
  def kt(u,v,up):
    f=k*(u-up)
    if (-fy)<f<fy:               return k,up
    if fy<=f    and v>0:  up=u-uy;return 0,up
    if fy<=f    and v<=0: up=u-uy;return k,up
    if f<=(-fy) and v<0:  up=u+uy;return 0,up
    else:                 up=u+uy;return k,up
  return kt

# Exercise #3 from 1st homework 2012
# 
mass =  1200.00       # kilograms
k  =   50000.00       # Newtons per metre
zeta =     0.05       # zeta is the damping ratio
fy =    3200.00       # yelding force, Newtons
t0 =       0.00
t1 =       0.50       # 
p0 =    1900.00       # 
p1 =    2100.00
uy =    fy/k          # yelding displacement, metres

# using the above constants, define the loading function
p=make_p(p0,p1,t0,t1)

# the following function, given the final displacement, the final
# velocity and the initial plastic deformation returns a) the tangent
# stiffness b) the final plastic deformation
kt=make_kt(k,fy) 

# we need the damping coefficient "c", to compute its value from the
# damping ratio we must first compute the undamped natural frequency 
wn=sqrt(k/mass)       # natural frequency of the undamped system
damp=2*mass*wn*zeta         # the damping coefficient

# the time step
h=0.02
# required duration for the response
t_end = 2.0
# the number of time steps to arrive at t_end
nsteps=int((t_end+h/100)/h)+1
# the maximum number of iterations in the Newton-Raphson procedure
maxiters = 30
# using the constant acceleration algorithm
# below we define the relevant algorithmic constants
gamma=0.5
beta=1./4.
gb=gamma/beta
a=mass/(beta*h)+damp*gb
b=0.5*mass/beta+h*damp*(0.5*gb-1.0)

# initialization of the system state and a bit more
t0=0.0
u0=0.0
up=0.0
v0=0.0
p0=p(t0)
(k0, up)=kt(u0,v0,up)
a0=(p0-damp*v0-k0*(u0-up))/mass

# print an header for the results
print "# Time     Load |     Disp.     Vel.   Accel.  |  Tang.St. Pl. Def. El.Force"

for i in range(nsteps):
    # output values at the beginning of each step
    print fmt % (t0,p0,u0,v0,a0,k0,up,k*(u0-up))
    # advance time, next external load value, etc
    t1 = t0 + h
    p1 = p(t1)
    Dp = p1 - p0
    Dp_= Dp + a*v0 + b*a0
    k_ = k0 + gb*damp/h + mass/(beta*h*h)
    # we prepare the machinery for the modified Newton-Raphson
    # algorithm.  if we have no state change in the time step, then the
    # N-R algorithm is equivalent to the standard procedure
    u_init=u0; v_init=v0 # initial state
    f_spring=k*(u0-up)   # the force in the spring
    DR=Dp_               # the unbalanced force, initially equal to the
                         # external load increment
    for j in range(maxiters):
        Du=DR/k_           # the disp increment according to the initial stiffness
        u_next = u_init + Du
        v_next = v_init + gb*Du/h - gb*v_init + h*(1.0-0.5*gb)*a0
                           # we are interested in the total plastic elongation
        oops,up=kt(u_next,v_next,up)
                           # because we need the spring force at the end
                           # of the time step
        f_spring_next=k*(u_next-up)
                           # so that we can compute the fraction of the
                           # incremental force that's equilibrated at the
                           # end of the time step
        df=f_spring_next-f_spring+(k_-k0)*Du
                           # and finally the incremental forces unbalanced
                           # at the end of the time step
        DR=DR-df
                           # finish updating the system state
        u_init=u_next; v_init=v_next; f_spring=f_spring_next
                           # if the unbalanced load is small enough (the
                           # criteria used in practical programs are
                           # energy based) exit the loop - note that we
                           # have no plasticization/unloading DR==0 at the
                           # end of the first iteration
        if abs(DR)<fy*1E-6 :
            break
        else:
            print "#", i, j, u_next, DR/fy
    # now the load increment is balanced by the spring force and
    # increments in inertial and damping forces, we need to compute the
    # full state at the end of the time step, and to change all
    # denominations to reflect the fact that we are starting a new time step.
    Du=u_init-u0
    Dv=gamma*Du/(beta*h)-gamma*v0/beta+h*(1.0-0.5*gamma/beta)*a0
    u1=u0+Du ; v1=v0+Dv
    k1,up=kt(u1,v1,up)
    a1=(p(t1)-damp*v1-k*(u1-up))/mass
    t0=t1; v0=v1; u0=u1 ; a0=a1 ; k0=k1 ; p0=p1

