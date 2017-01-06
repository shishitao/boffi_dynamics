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
Mn = [quad(lambda x,n=n: fi[n](x)*fi[n](x), 0., 1.)[0] for n in to8]
Gn = [Ln[n]/Mn[n] for n in to8] 

Vn = [Ln[n]*Gn[n] for n in to8] 
Hn = [quad(lambda x,n=n: Gn[n]*x*fi[n](x), 0., 1.)[0]/Vn[n] for n in to8]

print " n           &L_n         &M_n         &G_n         &V_n         &H_n        &BM_n\\\\"
for i in to8:
    print  "%2d  &%12f&%12f&%12f&%12f&%12f&%12f\\\\" % (i+1,Ln[i],Mn[i],Gn[i],Vn[i],Hn[i],Vn[i]*Hn[i])
