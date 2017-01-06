### as usual with python, import what we need need in our script 
from pylab import *
# <demo> stop
# we have a Tp=20s, and this period is divided in N=50 intervals
# the frequency interval is dw=2*pi/Tp=0.1256637 rad/s
# the time interval is dt=Tp/N=0.4s
# the Nyquist frequency is wny=2*pi*(1/Tp)*(N/2)=dw*(N/2)=dw*25
Tp = 20.0
N  = 50
dw = 2*pi/Tp
step = Tp/N
wny  = dw*25
# <demo> stop
# we need also a much higher sampling frequency, to simulate continuous functions
M=2000
# <demo> stop
# we want to compute the same functions on the same interval
# but with two different sampling rates
t_n=linspace(0.0,Tp,N+1)
t_m=linspace(0.0,Tp,M+1)
# <demo> stop
# we want to demonstrate by example the aliasing phenomenon, that is
# cos(n*w1*t)=cos((n-N)*t)
n = 31 ; print "n =",n, "     n-N =",n-N 
# the sampled functions are computed vector-wise
c1=cos(n*dw*t_n)
c2=cos((n-N)*dw*t_n)
# <demo> stop
# we want to plot also the continuous functions, so here we compute
# them vector-wise
c3=cos(n*dw*t_m)
c4=cos((n-N)*dw*t_m)
# and, at this point, we are ready to do our plots
# <demo> stop
figure(1) ; grid()
title(r'$\cos('+str(n)+'\,\omega_1t)$, continuous in red, 50 samples in blue')
figure(2) ; grid()
title(r'$\cos('+str(n-N)+'\,\omega_1 t)$, continuous in red, 50 samples in blue')
# <demo> stop
# we plot the (pseudo)continuous functions cos(n*dw*t) and cos((n-N)*dw*t
figure(1) ; plot(t_m,c3,'-r')
figure(2) ; plot(t_m,c4,'-r')
# <demo> stop
# the functions look different (one has visibly a higher frequency)
# but let's see what happens if we mark 50 equispaced points
# on both functions
figure(1) ; plot(t_n,c1,'ob')
figure(2) ; plot(t_n,c2,'ob')
# <demo> stop
# ok, the blue dots follow the same pattern in the two plots,
# it's time to plot the blue dots only, shifting n-N dots a little bit
figure(3) ; grid()
title('The two cosines, sampled at 50 points')
figure(3) ; plot(t_n,c1,'ob',t_n,c2-0.002*c2,'xr')
# <demo> stop
# we can zoom on a detail
axis([9.5,10.5,0.7,0.75])
# <demo> stop
show()
