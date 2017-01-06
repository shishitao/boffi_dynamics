from scipy import mat, sin, zeros
K = mat('1 0 0;0 2 0;0 0 3')
M = mat('4 1 0;1 4 1;0 1 2')/6.0
r = mat('0;-1;1')
def load(t): return r*sin(7.0*t) # 
h = 0.005; duration = 6.0
# linear acceleration coefficients
A = 3.0*M;   V = 6.0*M/h
Flex  = (K + 6.0*M/(h*h)).I
MI = M.I
# initial state
t = 0
x, v, p = mat(zeros((3,1))), mat(zeros((3,1))), r*sin(omega*t)
a = MI*(p - K*x)
print "%12.9f" % t, ' '.join(["%12.9f" % x_i for x_i in x])
# iterate with linear acceleration algorithm
while (t + h/2) < duration:
    t  = t + h
    dp = r*sin(omega*t) - p
    dx = Flex*(dp + A*a + V*v)
    dv = 3.0*dx/h - 3*v - a*h/2.0
    x, v, p = x + dx, v + dv, p + dp
    a  = MI*(p - K*x)
    print "%12.9f" % t, ' '.join(["%12.9f" % x_i for x_i in x])

