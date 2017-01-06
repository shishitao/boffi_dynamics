"""
This module computes the shock spectra for different types of
impulsive loadings, always under the assumption of zero damping.

Definitions:
t0      is the duration of the impulse,
T       is the natural period of vibration
ratio   is duration to period ratio
        ratio = t0 / T
beta is frequency ratio 
 beta = T  / (2 t0) ->  beta = 1 / 2a

alpha is an adimensional time,
   alpha = t/t0
 """
from math import sqrt, sin, asin, cos, pi

# rectangular is easy...

def rectangular_imp(ratio):
    """
    dummy docstring
    """
    if ratio < 0.5:
        return 2*sin(pi*ratio)
    else: return 2.0


def halfsine_imp(ratio):
    """Maximum response of a SDOF undamped system
    
This is harder, but for longer durations of the loading we have a
nice formula for the position of all extremes, so finding the
maximum is just a matter of appropriate bookkeping
"""
    
    def resp_ratio(alpha, ratio):
        """The response function in terms of
        alpha=t/t0, adimensional time parameter and
        ratio=t0/T,  duration/period ratio
        """
        dummy = 2*ratio*sin(alpha*pi)-sin(2*ratio*alpha*pi)
        dummy = 2*ratio*dummy
        return dummy/(4*ratio**2-1.)
        # return 2*ratio*(2*ratio*sin(alpha*pi)-
        # sin(2*ratio*alpha*pi) )/(4*ratio**2-1.)

    # short duration, maximum is "ballistic"
    if ratio < 0.5:
        return ratio*cos(ratio*pi) / (1./4 - ratio*ratio)
    # this must be special cased, as the prev. formula is
    # indeterminate for ratio=.5
    elif ratio == 0.5:
        return pi/2
    # duration is long enough to have a maximum during the forced response
    else:
        # we build a list of possible positions of stationary values,
        # starting with pos_1=1/(0.5+ratio)
        all_roots = [1/(0.5+ratio)]
        # for increasing values of n,
        for n in xrange(2, 100000):
            # until n is large with respect to ratio,
            # the duration/period ratio
            if 2*ratio <= n+n-1:
                break
            all_roots.append(-n/(0.5-ratio))
            all_roots.append(+n/(0.5+ratio))
        # when we have a list with the position of all stationary
        # values, we choose the maximum in absolute value and return
        # that value
        return max([resp_ratio(root, ratio) for root in all_roots], key=abs)



def triangular_imp(ratio):
    """From the particular integral that satisfies the EOM for an inverse ramp
                csi(t) = (P0/k)*(1-t/t0) = Dst * (1-t/t0)
we write (1) the general integral, (2) its time derivative, (3,4)
the initial conditions that give (5,6) the constants of integration
and finally (7) dividing by Dst we have the response function 
(1)                x(t)   = A sin(wn t) + B cos(wn t) + csi(t)
(2)                v(t)   = wn A cos(wn t) - B sin(wn t) - Dst/t0
(3)                x(0)   = B + Dst = 0
(4)                v(0)   = wn A - Dst/t0
(5)                     B = - Dst
(6)                     A = Dst / (wn t0)
(7)       R(t) = x(t)/Dst = sin(wn t)/(wn t0) - cos(wn t) + (1-t/t0)
with the usual substitutions, t=alpha*t0, ratio=t0/T, wn=2pi/T, we have
in (8) the Response function in terms of ratio, alpha, and in (9) its
time derivative
(8) R(alpha,ratio) = sin(2*ratio*pi*alpha) / (2*ratio*pi) - cos(2*ratio*pi*alpha) - alpha + 1
(9)         d R / d alpha = cos(...) + 2 ratio pi sin(...) - 1

by observation (i.e., plotting R in 0<alpha<1 for many values of ratio)
one knows that the greatest max is the first one,
so we have to find the smallest positive root of the derivative in 0<alpha<1

with w = 2 pi t0_t,

(10)        v(a) = 0 <=> cos(w a) + w sin(w a) = 1
                         sqrt(1+w^2) sin(w a - theta) = 1
with w b = w a - theta
                         sin(w b) = 1/sqrt(1+w^2)
                         w b = asin(1/sqrt(1+w^2)) = asin(k)

if, e.g, we take k=0.5, asin(0.5) is pi/6, but we are not interested
in the first root, we are interested in the second root,
                         w b_0 = pi-pi/6 [in our example]

from the plot it is apparent that theta = -asin(k), so 

                         w a_0 = w b_0 + theta = w b_0 - pi/6,


                +                         sin(w a) ******  
      1 ++ . . ... . . . . .*********. . . . . ... . . ++
        |       :       *               *       :       |
        |       :     *                   *     :       |
        |       :   *                       *   :       |
        |       : *                           * :       |
        |       *                               *       |
  k=0.5 ++ . . *:. . . . . . . . . . . . . . . .:* . . ++
        |    *  :                               :  *    |
        |   *   :                               :   *   |
        | *     :                               :     * |
        +*      +                               +      *+
      0 *-------+-------------------------------+-------*  tau = w b
        0     pi/6                            5pi/6    pi

                |-------------------------------|------------> w a
                0                             4pi/6
                                              w a_0
generalizing,

          w a_0 = pi - 2 asin(k)  ==>  a_0 = (pi - 2 asin(k))/w
    """
    if ratio < 0.371:
        # real impulsive response
        # 1. compute the displacement  and the velocity @ alpha=1
        disp = sin(2*ratio*pi) / (2*ratio*pi) - cos(2*ratio*pi)
        vel  = cos(2*ratio*pi) + 2*ratio*pi*sin(2*ratio*pi) - 1
        # return the srss
        return sqrt(disp*disp+vel*vel)
    else:
        omega = 2*pi*ratio
        root = (pi-2*(asin(1 /(omega*sqrt(1+1/omega/omega)))))/omega
        # and finally we return the value of R(alpha_0,ratio)
        return sin(omega*root)/omega - cos(omega*root) - root + 1
########################################################################
def main(command="print"):
    "Dummy doc string"
    n_iter = 1000
    if command != "print":
        import pylab as pl
        from pylab import vectorize as vcz
	from matplotlib import rcParams
        rcParams['figure.subplot.bottom'] = 0.15
        ratios = pl.linspace(5./n_iter,5,n_iter)
        pl.plot(ratios,vcz(rectangular_imp)(ratios))
        pl.plot(ratios,vcz(triangular_imp)(ratios))
        pl.plot(ratios,vcz(halfsine_imp)(ratios))
        pl.ylim(ymax=2.5)
        pl.grid(alpha=0.4)
        pl.xlabel(r"$\frac{t_o}{T_n}$",size=16)
        pl.xlabel(r"$t_o/T_n$",size=16)
        pl.ylabel("Peak Resp. Ratio")
        pl.xticks([ 0,   .371,   .5,   1,  2.5,   4.5,   5],
                  ["0","0.371","0.50","1","2.50","4.50","5"],
                  rotation=90)
        pl.show()
    else:
    # prints peak response for t_0/T_n ratios up to 5.0
    # skips t_0/T_n=0 because its leads to a division by zero
        for i in range(n_iter):
            ratio = 5./n_iter*(i+1)
            print "%10.6f     %12.10f  %12.10f  %12.10f" % (
                ratio,
                rectangular_imp(ratio),
                halfsine_imp(ratio),
                triangular_imp(ratio))

if __name__ == "__main__":
    main("plot")
