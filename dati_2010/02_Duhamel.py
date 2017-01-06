"""Computes the response of an undamped oscillator to a triangular load.

It is written with a bit of generality, so that it can be extended to
different type of loadings with little effort.
"""

from math import sin, cos, pi

def linear(y_0, y_1, x_0, x_1, k):
    """Returns the particular integral and its derivative for linear load.
    
    y_0, y_1 the ordinates @ the extremes of the loading,
    x_0, x_1 the abscissae @ the extremes of the loading,
    k      the stiffness of the oscillator.

    The returned functions are defined in terms of a new variable,
        csi=x-x_0,
    in the interval 0 <= csi <= x_0.
    """
    y_0, y_1 = y_0/k, y_1/k
    return (lambda x: y_0 + (y_1-y_0)/(x_1-x_0)*x,
            lambda x: (y_1-y_0)/(x_1-x_0))

def disp_vel(x_0, v_0, p_x, p_v, w_n):
    """Returns the displacement and velocity of the oscillator.

    x_0, v_0 the initial conditions,
    px, pv the particular integral and its time derivative,
    wn     the natural frequency of the oscillator.
    """
    A = (v_0 - p_v(0))/w_n
    B = x_0 - p_x(0)
    return (lambda t: A*sin(w_n*t)+B*cos(w_n*t)+p_x(t),
            lambda t: w_n*A*cos(w_n*t)-w_n*B*sin(w_n*t)+p_v(t))

# definition of the undamped oscillator, in terms of its mass and its
# natural period

MASS    =  600000.0
PERIOD  =       0.6

# derived constants
NAT_FREQ    = 2*pi/PERIOD
STIFFNESS   = NAT_FREQ**2*MASS

# definition of loading
TIME       = [  0.,      1.,  3., 5.]
FORCE      = [  0., 400000.,  0., 0.]
STEPS      = [ 100,     200, 200, 0]
PARTICULAR = [ linear, linear, linear, linear]

# initial conditions
disp_0 = 0.0
vel__0 = 0.0

for index in range(len(TIME)-1): # iterate over TIME intervals
    t_0, t_1  =  TIME[index],  TIME[index+1]
    p_0, p_1  = FORCE[index], FORCE[index+1]
    csi_func  = PARTICULAR[index]
    n_steps   = STEPS[index]
    px, pv    = csi_func(p_0, p_1, t_0, t_1, STIFFNESS)
    disp, vel = disp_vel(disp_0, vel__0, px, pv, NAT_FREQ)
    time_step = (t_1-t_0)/n_steps
    for time_index in range(n_steps): # iterate inside time intervals
        now = time_step*time_index
        print "%+12.4f %+12.7f %+12.7f" % (now+t_0, disp(now), vel(now))
    now = t_1-t_0
    disp_0, vel__0 = disp(now), vel(now)

print "%+12.4f %+12.7f %+12.7f" % (now+t_0, disp_0, vel__0)
