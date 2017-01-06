from math import *

def resp_elas(m,c,k, cC,cS,w, F, x0,v0):

    wn2=k/m ; wn=sqrt(wn2) ; beta=w/wn
    z=c/(2*m*wn)
    wd=wn*sqrt(1-z*z)

    det=(1.-beta**2)**2+(2*beta*z)**2
    R=((1-beta**2)*cS + (2*beta*z)*cC)/det/k
    S=((1-beta**2)*cC - (2*beta*z)*cS)/det/k
    D=F/k

    A=x0-S-D
    B=(v0+z*wn*A-w*R)/wd

    def x(t):
        return exp(-z*wn*t)*(A*cos(wd*t)+B*sin(wd*t))+R*sin(w*t)+S*cos(w*t)+D

    def v(t):
        return (-z*wn*exp(-z*wn*t)*(A*cos(wd*t)+B*sin(wd*t))
                  +wd*exp(-z*wn*t)*(B*cos(wd*t)-A*sin(wd*t))
                  +w*(R*cos(w*t)-S*sin(w*t)))

    return x,v

def resp_yeld(m,c, cC,cS,w, F, x0,v0):

    alpha=F/c
    det=w**2*(c**2+w**2*m**2)
    R=(+w*c*cC-w*w*m*cS)/det
    S=(-w*c*cS-w*w*m*cC)/det

    A=m*(w*R+alpha-v0)/c
    B=x0-A-S

    def x(t):
        return A*exp(-c*t/m)+B+R*sin(w*t)+S*cos(w*t)+alpha*t

    def v(t):
        return -c*A*exp(-c*t/m)/m+w*R*cos(w*t)-w*S*sin(w*t)+alpha

    return x,v

def bisect(f,t,x0,x1,count):
    h=(x0+x1)/2.0
    fh=f(h)-t
    if abs(fh)<1e-8 : return h
    f0=f(x0)-t
    if f0*fh>0 : return bisect(f,t,h,x1,count+1)
    return bisect(f,t,x0,h,count+1)

# from exercise 5.5 in Chopra's "Dynamics of Structures"
# data in imperial units, base units are kips, seconds and inches
mass=0.253302959106
k=10.0
zeta=0.05          
wn=sqrt(k/mass)    
damp=2*mass*wn*zeta
xy=0.75            
fy=k*xy            

t1=0.6
w=pi/t1
Po=10.0

x0=0.0
v0=0.0

x_next,v_next=resp_elas(mass,damp,k, 0.0,Po,w, 0.0, x0,v0)

t_yeld=bisect(x_next,xy,0.0,t1,1)

dtau=(t_yeld-0.0)/100.
for i in range (0,101):
    tau=dtau*i
    t=tau+0.0
    print t, x_next(tau), v_next(tau)
print "# yielding @ t =", t_yeld

x0=x_next(tau)
v0=v_next(tau)
cS=cos(w*tau)*Po
cC=sin(w*tau)*Po
x_next,v_next=resp_yeld(mass,damp,    cC,cS,w, -fy, x0,v0)

dt=(t1-t_yeld)/100.
for i in range(101):
    tau=i*dt
    t=tau+t_yeld
    print t, x_next(tau), v_next(tau)
print "# zero load @ t =", t1

cS=0.0 ; cC=0.0

x0=x_next(tau)
v0=v_next(tau)
x_next,v_next=resp_yeld(mass,damp,    cC,cS,w, -fy, x0,v0)

t2=t1+bisect( v_next, 0.0, 0, 0.3, 1.0)

dt=(t2-t1)/100
for i in range(101):
    tau=i*dt
    t=tau+t1
    print t, x_next(tau), v_next(tau)
print "# v=0, unloads @ t =", t2

x0=x_next(tau) ; v0=0.0
x_next,v_next=resp_elas(mass,damp,k, 0.0,0.0,w, k*x0-fy, x0,v0)

dt=(2.0-t2)/500
for i in range(501):
    tau=i*dt
    t=tau+t2
    print t, x_next(tau), v_next(tau)
