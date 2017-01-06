from scipy import *
from scipy.optimize import newton
from scipy.integrate import quad

to8=range(8)

betas=[newton(lambda x: cos(x)+1/cosh(x), (2*n+1)*pi/2) for n in to8]

coeff=[(cosh(b)+cos(b))/(sinh(b)+sin(b)) for b in betas]

fi=[lambda x, n = n:
    cosh(betas[n]*x)-cos(betas[n]*x) - coeff[n]*(sinh(betas[n]*x)-sin(betas[n]*x))
    for n in to8]

fi1=[fi[n](1.0) for n in to8]

fi=[(lambda x, n = n: (
    cosh(betas[n]*x) - cos(betas[n]*x) - coeff[n]*(sinh(betas[n]*x)-sin(betas[n]*x)))/fi1[n]
     ) for n in to8]

Ln = [quad(fi[n], 0., 1.)[0] for n in to8]
Mn = [quad(lambda x: fi[n](x)*fi[n](x), 0., 1.)[0] for n in to8]
Gn = [Ln[n]/Mn[n] for n in to8] 

print "    x    "+"G=%+5.3f "*8 % tuple(Gn)
for i in range (101):
    x=0.01*i
    print "%4.3f   " % (x,),
    for j in to8:
        print "%8.5f" % (0.8*Gn[j]*fi[j](x)+j+1,),
    print
