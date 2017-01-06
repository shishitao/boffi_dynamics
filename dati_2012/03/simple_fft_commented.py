"""Naive implementations of direct and inverse FFT

Uses lists instead of arrays, no optimizations
other than precomputing twiddles, no sanity tests
on inputs, just plain vanilla Decimation in Time..."""

from cmath import exp, pi

def d_fft(x,n):
    """Computes the direct fft of x, a list of n (n=2**k) complex values"""
    return _fft(x,n,[exp(-2*pi*1j*k/n) for k in range(n/2)])

def i_fft(x,n):
    """Computes the inverse fft of x, a list of n (n=2**k) complex values"""
    return [x/n for x in _fft(x,n,[exp(2*pi*1j*k/n) for k in range(n/2)])]

def _fft(x_vec, n_x, twiddle):
    """Implementation of Decimation in Time FFT, to be called by d_fft and i_fft.
    
    x_vec   is the signal to transform, a list of complex values
    n_x     is its length, results are undefined if n_x is not a power of 2
    twiddle is a list of twiddle factors, to be used in the computation
           (twiddles are precomputed by the caller functions)

    returns a list of complex values, to be normalized in case of an inverse transform"""


    if n_x == 1: return x_vec # bottom reached, DFT of a length 1 vec x is x

    # call fft with the even and the odd coefficients in x
    # the results are the so called even and odd  DFT's
    y_0 = _fft(x_vec[0::2], n_x/2, twiddle[::2])
    y_1 = _fft(x_vec[1::2], n_x/2, twiddle[::2])

    # assemble the partial results "in_place":
    # 1st half of full DFT is put in even DFT, 2nd half in odd DFT
    for k in range(n_x/2):
        y_0[k], y_1[k] = y_0[k]+twiddle[k]*y_1[k], y_0 [k]-twiddle[k]*y_1[k]

    # concatenate the two halves of the DFT and return to caller
    return y_0+y_1

def main():
    """Run some test cases"""
    from cmath import cos, sin, pi
    
    def testit(title, seq):
        """utility to format and print a vector and the ifft of its fft"""
        l_seq = len(seq)
        print "-"*5, title, "-"*5
        print "\n".join(["%10.6f :: %10.6f, %10.6fj" % (a.real, t.real, t.imag)
          for (a, t) in zip(seq, i_fft(d_fft(seq, l_seq), l_seq))])

    length = 32

    testit("Square wave", [+1.0+0.0j]*(length/2) + [-1.0+0.0j]*(length/2))
    testit("Sine wave",   [sin((2*pi*k)/length) for k in range(length)])
    testit("Cosine wave", [cos((2*pi*k)/length) for k in range(length)])
    
if __name__ == "__main__":
    main()
