"dummy doc string"
from math import sin, cos, pi, sqrt, atan2

def main():
    """The response of an undamped oscillator to a piecewise linear load.
    
    It is written with a bit of generality, so that it can be extended to
    different type of loadings with little effort.
    """
    
    # An undamped oscillator is defined in terms of two parameters, here
    # the mass and the period, from which it is possible to derive the
    # other relevant parameters
    mass        = 600000.0
    period      = 0.6
    omega_n     = 2*pi/period
    stiffness   = omega_n**2*mass
    
    # Definition of loading: we have a piecewise, continuos linear
    # loading, so for each interval where the loading is described by the
    # same linear function we define the relevant parameters:
    times       = [0., 1., 3., 5.]
    forces      = [0., 400000., 0., 0.] # Newtons
    n_steps     = [100, 200, 200]       # four time points, 3 intervals 
    # which function is used, in each interval, to get the particular solution
    force_types = [ linear, linear, linear] 
    # of course, if we want to use different types of loading, then we
    # have to modify the data structures, in particular the FORCE vector
    # must contain, for each interval, a properly conceived vector that
    # describes the (more complex) loading
    
    # initial conditions
    disp_0 = 0.0
    vel__0 = 0.0
    
    print times
    for index in range(len(times)-1):
        # unfold the relevant parameters from their containers
        t_0, t_1  =  times[index],  times[index+1]
        p_0, p_1  = forces[index], forces[index+1]
        part_sol  = force_types[index]
        steps     = n_steps[index]
        time_step = (t_1-t_0)/steps
        # compute the FUNCTIONS that give the solution in the current interval
        # NB: these are functions of tau=t-t_0
        xi_dis,   xi_vel  = part_sol(p_0, p_1, t_0, t_1, stiffness)
        disp, vel = disp_vel(disp_0, vel__0, xi_dis, xi_vel, omega_n)
        # do a loop printing the values of t=tau+t_0, disp, vel
        for time_index in range(steps): # NB time_index starts from zero
            tau = time_step*time_index
            print "%+12.4f %+12.7f %+12.7f" % (tau+t_0, disp(tau), vel(tau))
        # now t=t_1-time_step
        # we need to compute the disp and vel at tau=t_1-t_0 to set
        # the initial values for the following interval
        tau = t_1-t_0
        disp_0, vel__0 = disp(tau), vel(tau)
    # we have done all intervals, but we have not printed yet the last
    # time point... let's do it
    print "%+12.4f %+12.7f %+12.7f" % (tau+t_0, disp_0, vel__0)
# END OF MAIN

def linear(f_0, f_1, t_0, t_1, k):
    r"""Returns \xi(t), \dot\xi(t) for an undamped SDOF subject to linear load.
    
    f_0, f_1 the initial and final values of the linear force 
    t_0, t_1 the initial and final values of the time
    k      the stiffness of the oscillator.

    The returned functions are defined in terms of a new variable,
                         tau=t-t_0,
    in the interval 0 <= tau <= t_1-t_0.
    """
    x_0, x_1 = f_0/k, f_1/k
    delta_x  = x_1 - x_0
    delta_t  = t_1 - t_0
    velocity = delta_x / delta_t
    return (lambda tau: x_0 + velocity*tau,
            lambda tau: velocity)
# END OF LINEAR

def disp_vel(x_0, v_0, p_x, p_v, w_n):
    """Returns the displacement and velocity of the oscillator.

    x_0, v_0 the initial conditions,
    px, pv the particular integral and its time derivative,
    wn     the natural frequency of the oscillator.
    """

    A     = (v_0 - p_v(0))/w_n
    B     = x_0 - p_x(0)
    D     = sqrt(A*A + B*B)
    theta = atan2(-B, A)

    return (lambda t:     D*sin(w_n*t-theta)+p_x(t),
            lambda t: w_n*D*cos(w_n*t-theta)+p_v(t))
# END OF DISP_VEL

if __name__ == "__main__":
    main()
