from scipy import *
from scipy.linalg import det, eigh
from sys import exit
def ltxp(data,name,fmt="%10.4f",title="",pre="",post=""):
    delim={"mat":"b", "vet":"B", "det":"V", "norm":"v"}
    if title: print "% ----- "+title+" -----"
    print pre+"\\begin{"+delim[name]+"matrix}"
    print "\\\\\n".join(["&".join(map(lambda x: fmt%(x,),line)) for line in asarray(data)])
    print "\\end{"+delim[name]+"matrix}"+post

def integ(a, b, l):
    integral = polyint(a*b)
    return integral(l) -integral(0)

p = poly1d; pi = polyint

m = [[p([1, 0]), p([1, 1]), p([1, 2])],
     [p([0, 0]), p([0, 0]), p([1, 0])],
     [p([0, 0]), p([1, 0]), p([1, 1])]]
l = [1,1,1]

flex = matrix(zeros((3,3)))

for i in (0, 1, 2):
    for j in (0, 1, 2):
        flex[i,j] = sum([integ(a, b, lg) for a, b, lg in zip(m[i], m[j], l)])

K = flex.I
M = matrix("1. 0. 0. ; 0. 1. 0.; 0. 0. 2.")

ltxp(6*flex,"mat",
     title="The flexibility matrix normalized by (L^3)/(6*EJ)",
     pre=r"\bm{F}=\frac{L^3}{6EJ}")

ltxp(13*K/3, "mat",
     title="The stiffness matrix normalized by (3EJ)/(13L^3)*",
     pre=r"\bm{K}=\frac{3 EJ}{13 L^3}")


ltxp(M,"mat",
     title="The mass matrix: m*",pre=r"\bm{M}=m")

K_xx = K[:2,:2] ; K_xg = K[:2,2:] ; E = -K_xx.I*K_xg
ltxp(-16*(K[:2,:2]).I*K[:2,2:], "mat",
     title='The influence matrix "E"', pre=r'\bm{E}=\frac{1}{16}')
exit()

evals, evecs = eigh(K,M)
#print "The eigenvalues:\n", evals,","
#print "the eigenvectors:\n", evecs,"."

ltxp(evals.reshape(1,3),"vet",
     title="The eigenvalues",pre=r"\bm{\Omega}=",post="^{T}")
ltxp(evecs,"mat",title="The eigenvectors",pre=r"\bm{\Psi}=")
ltxp(evecs.T*M*matrix("0;0;1"),"mat",
     title="Modal load vector, to multiply \ddot{u}_g",
     pre=r"\bm{p}_\text{eff}=m",post=r"\ddot u_\text{g}")
