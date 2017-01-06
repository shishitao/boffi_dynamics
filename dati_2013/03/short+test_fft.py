from cmath import exp, cos, sin, pi 
def f(x,n,w): return (lambda y=f(x[::2],n/2,w[::2]),z=f(x[1::2],n/2,w[::2]):reduce(lambda x,y:x+y,zip(*[(y[k]+w[k]*z[k],y[k]-w[k]*z[k]) for k in range(n/2)])))() if n>1 else x

def dfft(x,n): return               f(x,n,[exp(-2*pi*1j*k/n) for k in range(n/2)]) 
def ifft(x,n): return [x/n for x in f(x,n,[exp(+2*pi*1j*k/n) for k in range(n/2)])] 

def main(): 
    """Run some test cases""" 
    def testit(title, seq): 
        """utility to format and print a vector and the ifft of its dfft""" 
        l_seq = len(seq) 
        print "-"*5, title, "-"*5 
        print "\n".join(["%10.6f :: %10.6f, %10.6fj" % (a.real, t.real, t.imag) 
          for (a, t) in zip(seq, ifft(dfft(seq, l_seq), l_seq))]) 
    length = 32 
    testit("Square wave", [+1.0+0.0j]*(length/2) + [-1.0+0.0j]*(length/2)) 
    testit("Sine wave",   [sin(2*pi*k/length) for k in range(length)]) 
    testit("Cosine wave", [cos(2*pi*k/length) for k in range(length)]) 
# ---------------------------------------------------------------------- 
if __name__ == "__main__": 
    main() 

