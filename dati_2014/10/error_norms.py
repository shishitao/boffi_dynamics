import scipy as sp
from scipy import mat, sqrt
from scipy.linalg import eigh

def lp(data,name,fmt="%+10.4f",title=""):
    delim={"mat":"b", "vet":"B", "det":"V", "norm":"v"}
    if title: print "% ----- "+title+" -----"
    print "\\begin{"+delim[name]+"matrix}"
    print "\\\\\n".join(["&".join(map(lambda x: fmt%(x,),line)) for line in sp.asarray(data)])
    print "\\end{"+delim[name]+"matrix}"
    
def derritz(n,r):
    global T, phi

    phi=mat(sp.zeros((n,n)))
    
    y=F*r
    b=sqrt((y.T*y)[0,0])
    y=y/b
    phi[:,0]=y[:,0]

    for i in range(1,n):
        j=i-1
        y=F*phi[:,j]
        alpha=y.T*phi
        for k in range(i):
            y=y-alpha[0,k]*phi[:,k]
        b=sqrt((y.T*y)[0,0])
        phi[:,i]=y[:,0]/b
    T=phi.T*F*phi
    T=sp.diagflat(sp.diag(T))+sp.diagflat(sp.diag(T,1),1)+sp.diagflat(sp.diag(T,-1),-1)

# Fixed data

n=5

M=mat(sp.eye(n))
K=mat(sp.eye(n))*2
K[n-1,n-1]=1.
for i in range(n-1):
    K[i+1,i]=-1.
    K[i,i+1]=-1.
F=K.I

lp(M,'mat',title="Mass matrix (=I) M")
lp(K,'mat',title="Stiffness matrix K")
lp(F,'mat',title="Flex.lity matrix F")

# evecs are normalized with respect to M
evals, evecs = eigh(K,M)
L=mat(sp.diagflat(evals))
lp(L,'mat',title="Eigenvalues Matrix L")
lp(evecs,'mat',title="Eigenvectors matrix, \Psi")

for r in ( sp.mat((0,0,0,0.,1.)).T,sp.mat((0,0,0,-2.,1.)).T, sp.mat((1,1,1,1.,1.)).T):
    lp(r.T,'vet',title="Load vector transposed r^T")
    derritz(n,r)
    lp(phi,'mat',title="Derived Ritz Vectors matrix \Phi")
    # lp(T,'mat',title="Tridiagonal DRV matrix, T")

    Gamma=evecs.T*r
    # lp(Gamma.T,'vet',title="Modal partecipation factors")
    Gammh=phi.T*r
    # lp(Gammh.T,'vet',title="DRV's partecipation factors")


    f_m=evecs*sp.diagflat(Gamma)
    f_r=phi*sp.diagflat(Gammh)

    # lp(f_m,'mat',title="Modal forces matrix")
    # lp(f_r,'mat',title="DRV's forces matrix")

    den=sp.dot(r.T,r)[0,0]
    e_m = r ; e_r = r
    for i in range(n):
        e_m = e_m-f_m[:,i] ; e_r = e_r-f_r[:,i] ;
        print "%3d   %10.7f   %10.7f" % (i+1, sp.dot(r.T,e_m)[0,0]/den,  sp.dot(r.T,e_r)[0,0]/den )
