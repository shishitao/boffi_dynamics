from scipy import *
from scipy.linalg import eigh

def ltxp(data,name,fmt="%10.4f",title="",pre="",post=""):
    """print a matrix in LaTeX format"""
    
    delim={"mat":"b", "vet":"B", "det":"V", "norm":"v"}
    if title: print "% ----- "+title+" -----"
    print pre+"\\begin{"+delim[name]+"matrix}"
    print "\\\\\n".join(["&".join(map(lambda x: fmt%(x,),l)) for l in asarray(data)])
    print "\\end{"+delim[name]+"matrix}"+post

def coef_by_coef(a,b):
    "Given 2 matrices gives the coef x coef product (not the matrix product)"
    return asarray(a).flatten()*asarray(b).flatten()

K = matrix('1 0 0;0 2 0;0 0 3')
M = matrix('4 1 0;1 4 1;0 1 2')
p = matrix('0;-1;1')

ltxp(K,'mat',title='stiffness',pre=r'\bm{K}=k')
ltxp(M,'mat',title='mass',pre=r'\bm{M}=\frac{m}{6}')
M = M/6.
ltxp(p.T,'vet',title='external load',pre=r'\bm{p}=\frac{kL}{200}',post='^T')

evals, evecs = map(matrix,eigh(K,M))
Lambda = matrix(diag(asarray(evals).flatten()))
ltxp(evals,'vet',title='normalized eigenvalues',pre=r'\bm\Lambda=',post='^T')
ltxp(evecs,'mat',fmt="%+12.8f",title='eigenvector matrix',pre=r'\bm\Psi=')

# s-s response, ss is normalized with respect to L/200
ss = (K-49*M).I*p

# initial velocities in modal coordinates, normalized as above
dotq0 = -evecs.T*M*7.0*ss

# divide the initial velocities by w_i to get the sine coefficients, normalized
b =  sqrt(Lambda.I)*dotq0

ltxp(ss,'mat',fmt="%+12.8f",title='s-s nodal responses',pre=r'\bm\xi=')
ltxp(dotq0,'vet',fmt="%+12.8f",title='modal vel. initial values')
ltxp(b,'mat',fmt="%+12.8f",title='sine coeff in modal coor.',pre=r'\bm\b=\frac{L}{200}')

# below a too smart way of writing 
x3_by_qs = matrix(coef_by_coef(evecs[2],b.T))
ltxp(x3_by_qs,'vet',fmt="%+12.7g",title='sine coefficients in structural coordinates')
ltxp(sqrt(evals),'vet',fmt="%+12.7g",title="frequencies")
