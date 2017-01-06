from __future__ import division
import sympy as sy
from sympy.printing import latex
# ----------------------------------------------------------------------
def place(i,j,k,l, K,a):
    "Ad hoc placement for our simple problem"
    K[i,i] += a[0,0]; K[i,j] += a[0,1]; K[i,k] += a[0,2]; K[i,l] += a[0,3]
    K[j,i] += a[1,0]; K[j,j] += a[1,1]; K[j,k] += a[1,2]; K[j,l] += a[1,3]
    K[k,i] += a[2,0]; K[k,j] += a[2,1]; K[k,k] += a[2,2]; K[k,l] += a[2,3]
    K[l,i] += a[3,0]; K[l,j] += a[3,1]; K[l,k] += a[3,2]; K[l,l] += a[3,3]
def fx(f):
    "Make a function that can be used by Matrix() to fill K_element"
    def inte(i,j):
        return sy.integrate(f[i]*f[j],(x,0,L))   
    return inte
def partition(mtx, n):
    "Returns K_aa, K_ab, K_ba, K_bb, with K_aa dimensions n by n"
    return mtx[:n, :n], mtx[:n, n:], mtx[n:, :n], mtx[n:, n:]
# ----------------------------------------------------------------------

x, m, L, EJ = sy.symbols('x m L EJ')

# beam finite element shape functions
f =  [    +1       -3*(x/L)**2 +2*(x/L)**3,
                   +3*(x/L)**2 -2*(x/L)**3,
       L*(+1*(x/L) -2*(x/L)**2 +1*(x/L)**3),
       L*(         -1*(x/L)**2 +1*(x/L)**3)]

# comp is a function that given i,j computes \int_0^L f"_i f"_j dx
comp=fx([sy.diff(f_x,x,2) for f_x in f])

# k_el is normalized with respect to EJ (implicitly) and L**3
k_el=L**3*sy.Matrix(4,4,comp)

# create and fill the 10x10 stiffness matrix, first 5 translational DOFs,
#   then 5 rotational DOFs
K=sy.zeros(10)
for n in range(4):
    place(n,n+1,n+5,n+6, K,k_el)

# partition K and do static condensation
Kdd, Kdr, Krd, Krr = partition(K, 5)
d2r = Krr.inv()*Krd # rotations = - d2r * displacements
Kred=Kdd-Kdr*d2r    # Kred is for both nodes and supports
Ksav=Kred[:,:]
# x_1 left supp., x_2 left node, x_3 centre supp., x_4 right node,
# x_5 right supp.
# rearrange stiffness so that nodes' DOFS come first
Kred.row_swap(0,1);Kred.row_swap(2,3);Kred.row_swap(1,2)
Kred.col_swap(0,1);Kred.col_swap(2,3);Kred.col_swap(1,2)
# x_1, x_2 left right nodes, then left, centre and right supports

# partition Kred to get s_tructural and g_round matrices
Kss, Ksg, Kgs, Kgg = partition(Kred, 2)
Emat=-Kss.inv()*Ksg

# eigenvectors problem is Kss psi = w^2/(EJ/mL^3) I psi
ee2,ee1=Kss.eigenvects()

Psi = ee1[2][0]
Psi = Psi.col_insert(1,ee2[2][0])
M_star=Psi.transpose()*Psi
Ell=Psi.transpose()*Emat
Gamma=M_star.inv()*Ell

print latex(Psi)
print latex(M_star)
print latex(M_star.inv()*Ell)
gg=Gamma.evalf()
print latex(gg)

D10,D11,D12,D20,D21,D22 = sy.symbols("D11 D12 D13 D21 D22 D23")
D=sy.Matrix(((D10, D11,D12),(D20,D21,D22))).transpose()
GammaD=sy.zeros((2,1))
GammaD[0,0]=(Gamma[0,:]*D[:,0])[0,0]         
GammaD[1,0]=(Gamma[1,:]*D[:,1])[0,0]
psi11, psi21, psi12, psi22 = sy.symbols("psi11 psi21 psi12 psi22")
psi=sy.Matrix((((psi11,psi12)),(psi21,psi22)))
print latex(psi*GammaD)





print r"\documentclass[12pt]{article}\usepackage{amsmath}\usepackage{bm}\begin{document}"
print (r"The stiffness matrix for the 10x10 model is\
\[\bm K_{10\times10}=\frac{EJ}{L^3}\text{"+latex(K,mat_delim="[")+r"}\]" +"\n"+
r"the mapping between translational and rotational degrees of freedom is given by\
"+r"\[\vec{\bm\phi}=\frac{1}{56L}\text{"+latex(56*L*d2r,mat_delim="[")+r"}\vec{\bm x}.\]")
print
print r"Following static condensation and reordering rows and columns,\
the partitioned stiffness matrices are\
\[\bm K = \frac{EJ}{28L^3}\text{"+latex(28*Kss,mat_delim="[")+r"},\;\
\bm K_\text{g}=\frac{EJ}{28L^3}\text{"+latex(28*Ksg,mat_delim="[")+r"}\
\bm K_\text{gg}=\frac{EJ}{28L^3}\text{"+latex(28*Kgg,mat_delim="[")+r"}.\]"
print r"\
The influence matrix is\
\[ {\bm E} = {\bm K}^{-1}{\bm K}_\text{g}=\frac{1}{32}\text{"+latex(32*Emat,mat_delim="[")+r"}.\]"
print
print r"The eigenvector matrix is\
\[{\bm \Psi}=\text{"+latex(Psi,mat_delim="[")+r"}\]\
the matrix of modal masses is\
\[ {\bm M}^\star={\bm\Psi}^T {\bm M} {\bm \Psi} = m\text{"+latex(M_star,mat_delim="[")+r"}\]\
the matrix of the non normalized modal partecipation coefficients is\
\[ {\bm L} = {\bm\Psi}^T {\bm M} {\bm E} = m\text{"+latex(Ell,mat_delim="[")+r"}\]\
and, finally, the matrix of modal partecipation factors,\
\[ {\bm \Gamma} = ({\bm M}^\star)^{-1} {\bm L} = \text{"+latex(Gamma,mat_delim="[")+r"}\]"
print
print r"Denoting with $D_{ij}=D_{ij}(t)$ the response function for mode $i$ due to\
ground excitation $\ddot{x}_{\text{g}j}$, the response can be written\
\[ {\bm x}-{\bm E}{\bm x}_\text{g} = \text{"+latex(psi*GammaD)+r"} = \text{"+latex(Psi*GammaD)+r"}.\]"
print latex(GammaD)
print r"\end{document}"
