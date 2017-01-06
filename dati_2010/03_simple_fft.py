"""Naive implementations of direct and inverse FFT

Uses lists instead of arrays, no optimizations
other than precomputing twiddles, no sanity tests
on inputs, just plain vanilla Decimation in Time..."""

def fft(x_vec, n_x, sign=-1, twiddle=None):
    """x_vec is the signal to transform
n_x   is its length
sign  is the sign in the exponential of imaginary argument
      = -1 for direct transform (default)
      = +1 for inverse transform
      
if this function is used for inverse transform (use ifft instead)
the result must be normalized, i.e., divided by n_x"""

    from cmath import exp, pi
    
    if not twiddle:
        twiddle = [exp(sign*2*pi*1j*k/n_x) for k in range(n_x/2)]

    if n_x == 1:    # bottom reached, DFT of a length 1 vec x is x
        return x_vec

    # call fft with the even and the odd coefficients in x
    # to compute the so called even and odd  DFT's
    y_0 = fft(x_vec[0::2], n_x/2, sign, twiddle[::2])
    y_1 = fft(x_vec[1::2], n_x/2, sign, twiddle[::2])

    # assemble the partial results "in_place":
    # 1st half of full DFT is put in even DFT, 2nd half in odd DFT
    for k in range(n_x/2):
        y_0[k], y_1[k] = y_0[k]+twiddle[k]*y_1[k], y_0 [k]-twiddle[k]*y_1[k]

    # concatenate the two halves of the DFT and return to caller
    return y_0+y_1

def ifft(x_vec, n_x):
    """Define the inverse FFT in terms of the direct FFT"""
    return [y_vec/n_x for y_vec in fft(x_vec, n_x, sign=1)]

def main():
    """Run some test cases"""
    from cmath import cos, sin, pi
    
    def testit(title, seq):
        """utility to format and print a vector and the ifft of its fft"""
        l_seq = len(seq)
        print "-"*5, title, "-"*5
        print "\n".join(["%10.6f :: %10.6f, %10.6fj" % (a.real, t.real, t.imag)
          for (a, t) in zip(seq, ifft(fft(seq, l_seq), l_seq))])

    length = 32

    testit("Square wave", [+1.0+0.0j]*(length/2) + [-1.0+0.0j]*(length/2))
    testit("Sine wave",   [sin((2*pi*k)/length) for k in range(length)])
    testit("Cosine wave", [cos((2*pi*k)/length) for k in range(length)])
    
if __name__ == "__main__":
    main()
