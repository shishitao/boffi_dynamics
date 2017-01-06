from scipy import *
from scipy.linalg import eigh

story_stiffness = range(23,11,-1)
story_stiffness.append(0)

K = matrix(zeros((13,13)))
M = matrix(zeros((12,12)))

for i in range(12):
    M[i,i] = 1.0
    K[i,i] = story_stiffness[i] + story_stiffness[i+1]
    K[i,i+1] = -story_stiffness[i+1]
    K[i+1,i] = -story_stiffness[i+1]

K = K[:12,:12]

print "story stifnesses:"
print story_stiffness
print 
print "normalized mass matrix M/m:"
print M
print
print "normalized stiffness matrix K/k:"
print K
print

evals, evecs = eigh(K,M,eigvals=(0,3))
print "first four eigenvalues of the system"
print evals
print

phi=matrix("""
 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2;
 0.1, 0.2, 0.3, 0.4, 0.4, 0.3, 0.2, 0.1, 0.0,-0.1,-0.2,-0.3;
 0.1, 0.3, 0.4, 0.3, 0.1,-0.1,-0.3,-0.3,-0.1, 0.1, 0.2, 0.4;
 0.1, 0.3, 0.1,-0.1,-0.3,-0.3,-0.1, 0.1, 0.3, 0.1,-0.1,-0.4""").T

print "initial subspace vectors"
print phi
print

def iterate(phi):
    Mr = phi.T*M*phi
    Kr = phi.T*K*phi
    zval, zvec = eigh(Kr,Mr)
    return zval, phi*zvec, K.I*M*phi*zvec
