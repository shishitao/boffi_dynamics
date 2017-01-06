######################################################################
# Preliminaries,
import scipy as sp
mat=sp.matrix
from scipy.linalg import inv
######################################################################
# an utility function to format a matrix for inclusion in a LaTeX file
def latex_print(data,name,fmt="%10.4f",title=""):

    delim={"mat":"b",
           "vet":"B",
           "det":"V",
           "norm":"v"}

    if title:
        print "% ----- "+title+" -----"
    print "\\begin{"+delim[name]+"matrix}"
    print "\\\\\n".join(["&".join(map(lambda x: fmt%(x,),line)) for line in sp.asarray(data)])
    print "\\end{"+delim[name]+"matrix}"
######################################################################
Mass=mat(((2,0,0,),
          (0,3,0,),
          (0,0,4,),));
Mass=100000.*Mass
latex_print(Mass,"mat",title="Mass Matrix",fmt="%10.0f")
    
######################################################################
Stif=mat(((+1,-1,+0),
          (-1,+3,-2),
          (+0,-2,+5)))
Stif=120e6*Stif
latex_print(Stif,"mat",title="Stiffness Matrix",fmt="%10.0f")

######################################################################
#                    roots finds the roots of the poly defined by
#                          the list of coefficients (1, -11.4,m ...)
Omegas=mat(sorted(sp.roots((1,-11/4.,15/8.,-1/4.))))*1200.

Eigenv=mat(sp.zeros((3,3)))
# This sets the 0 row of the eigenv matrix to ones
Eigenv[0,:]=1.,1.,1.
# this is a {0,1} column vector
known=mat(((1,),(0,)))

# solve the eq. of free vibrations for psi_0i = 1
for i in range(3):
    Omega2=Omegas[0,i]/1200
    coef=mat(((3.-3.*Omega2,-2.),(-2.,5.-4.*Omega2)))
    bottom=coef.I*known
    # this sets the bottom part of each eigenvector
    Eigenv[1:,i]=bottom
latex_print(Eigenv,"mat",title="Eigenvectors Matrix")

MStar=Eigenv.T*Mass*Eigenv
latex_print(MStar,"mat",title="Modal Masses Matrix")

KStar=Eigenv.T*Stif*Eigenv
latex_print(KStar,"mat","%10.5e",title="Modal Stiffnesses Matrix")

MStar=Eigenv.T*Mass*Eigenv
latex_print(MStar/1000.,"mat",title="Modal Masses Matrix, in tons")

KStar=Eigenv.T*Stif*Eigenv
latex_print(KStar/1E6,"mat","%10.2f",title="Modal Stiffnesses Matrix, in MN/m")

q_0=MStar.I*Eigenv.T*Mass*mat((5,4,3)).T
latex_print(sp.mat(((5,4,3),)).T,"vet",title="Initial displacements, nodal coo.")
latex_print(q_0,"vet",title="Initial displacements, modal coo.")

qdot_0=MStar.I*Eigenv.T*Mass*mat((0,9,0)).T
latex_print(mat((0,9,0)).T,"vet",title="Initial velocities, nodal coo.")
latex_print(qdot_0,"vet",title="Initial velocities, modal coo.")

# q_i    =     A_i sin(w_i t) + B_i cos(w_i t)
# qdot_i = w_i(A_i cos(w_i t) - B_i sin(w_i t))

Bs=q_0
As=mat(sp.diagonal(qdot_0/sp.sqrt(Omegas))).T

latex_print(As,"vet",title="Sine coefficients for modal disp.s")
latex_print(Bs,"vet",title="Cosine coefficients for modal disp.s")


ampli=sp.real(sp.sqrt(sp.power(As,2)+(sp.power(Bs,2))))
phase=sp.arctan2(As,Bs)
latex_print(ampli,"vet",title="Cosine only amplitudes for modal disp.s")
latex_print(phase,"vet",title="Cosine only phases for modal disp.s")


# q_i(t) = ampli_i*cos(w_i-phase)

print "% Nodal displacements, in mm\n\\begin{align*}"
for i in range(3):
    print r"  x_%d & = " % (i+1,),
    for j in range(3):
        print r"%+6.3f \cos(%10.3f t %+10.3f) " % (Eigenv[i,j]*ampli[j], sp.sqrt(Omegas[0,j]), phase[j]),
    print r"\\"
print "\\end{align*}"

print "% Nodal forces, in kN\n\\begin{align*}"

for i in range(3):
    print r"x_%d & = " % (i+1,),
    for j in range(3):
        print r"%+6.3f \cos(%10.3f t %+10.3f) " % (Mass[i,i]*Omegas[0,j]*Eigenv[i,j]*ampli[j]/1E6, sp.sqrt(Omegas[0,j]), phase[j]),
    print r"\\"
print "\\end{align*}"

## half-sine
#t1=0.02 # seconds
#p=mat((2.5e6,5e6,5e6)).T # Newtons
## modal loads, normalized
#pl=MStar.I*Eigenv.T*p
##the impulse, and the final velocity, as pl was normalized, is
#qdot_0 = pl*t1/(sp.pi/2)
#print qdot_0, sp.diagonal(qdot_0/sp.sqrt(Omegas))
