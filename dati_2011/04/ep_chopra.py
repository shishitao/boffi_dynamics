"""Using data from an example in Chopra's "Dynamics of Structures"
this program integrates the equation of motion of an elastic-
perfectly plastic SDOF system, using the modified Newton-Raphson
method to achieve equilibrium during changes of stater."""

from math import pi, sqrt, sin

HEAD = "# Time     Load |     Disp.     Vel.   Accel.  "
HEAD += "|  Tang.St. Pl. Def. El.Force"
FMT = "%6.3f %8.4f    %+8.4f %+8.4f %+8.4f     %+8.4f %+8.4f %+8.4f"

def make_p(load_max, duration):
    """make_p(load_max,duration)
returns a 1/2 sine impulse load function, half_sine(t)"""
    def half_sine(time):
        "half-sine"
        if time < duration:
            return load_max*sin(time*pi/duration)
        else:
            return 0.0
    return half_sine

def make_kt(stiffness, yielding_force, yielding_disp):
    """make_kt(stiffness,yielding_force)
returns a function k_t(u,v,up) returning k_tan, up"""
    def k_t(disp, vel, plastic_def):
        "function kt_and_up(u, v, u_p) returning kt_and_up, u_p"
        force = stiffness*(disp-plastic_def)
        if (-yielding_force) < force < yielding_force:
            return stiffness, plastic_def
        if yielding_force <= force   and vel > 0:
            plastic_def = disp-yielding_disp
            return 0, plastic_def
        if yielding_force <= force   and vel <= 0:
            plastic_def = disp-yielding_disp
            return stiffness, plastic_def
        if force < (-yielding_force) and vel < 0:
            plastic_def = disp+yielding_disp
            return 0, plastic_def
        else:
            plastic_def = disp+yielding_disp
            return stiffness, plastic_def
    return k_t

# from exercise 5.5 in Chopra's "Dynamics of Structures"
# data in imperial units, base units are kips, seconds and inches

MASS = 0.253302959106
STIFFNESS = 10.0
# zeta is the damping ratio
ZETA = 0.05
# yelding displacement
U_YIELD = 0.75
# yelding force
F_YIELD = STIFFNESS*U_YIELD

# parameters of the loading function

# half-sine impulse duration
DURATION = 0.6
# half-sine impulse intensity
LOAD_MAX = 10.0

# we need the damping coefficient "DAMPING", to compute its value from the
# damping ratio we must first compute the undamped natural frequency 

# natural frequency of the undamped system
W_N = sqrt(STIFFNESS/MASS)
# the damping coefficient
DAMPING = 2*MASS*W_N*ZETA

# the time step used in Chopra's exercise was 0.1, a bit too long, 
STEP = 0.01
# the number of time steps to arrive at t = 1.0
NSTEPS = int((1.0+STEP/100)/STEP)+1
# the maximum number of iterations in the Newton-Raphson procedure
MAXITERS = 30

# the exercise in Chopra's uses the constant acceleration algorithm
# below we define the relevant algorithmic constants
GAMMA = 0.5
BETA = 1./4.
G_B = GAMMA/BETA
A = MASS/(BETA*STEP)+DAMPING*G_B
B = 0.5*MASS/BETA+STEP*DAMPING*(0.5*G_B-1.0)

def main():
    """integrate the equation of motion
    """
    # initialization of the system state and a bit more

    load = make_p(LOAD_MAX, DURATION)
    kt_and_up = make_kt(STIFFNESS, F_YIELD, U_YIELD) 

    t_0 = 0.0
    u_0 = 0.0
    u_p = 0.0
    v_0 = 0.0
    p_0 = load(t_0)
    k_0, u_p = kt_and_up(u_0, v_0, u_p)
    a_0 = (p_0-DAMPING*v_0-k_0*(u_0-u_p))/MASS

    # print an header for the results
    print HEAD


    for i in range(NSTEPS):
        # output values at the beginning of each step
        print FMT % (t_0, p_0, u_0, v_0, a_0, k_0, u_p, STIFFNESS*(u_0-u_p))
        # advance time, next external load value, etc
        t_1 = t_0 + STEP
        p_1 = load(t_1)
        delta_p = p_1 - p_0
        delta_p_eff =  delta_p + A*v_0 + B*a_0
        k_eff = k_0 + G_B*DAMPING/STEP + MASS/(BETA*STEP*STEP)
        # we prepare the machinery for the modified Newton-Raphson
        # algorithm.  if we have no state change in the time step, then
        # the N-R algorithm is equivalent to the standard procedure

        # initial state
        u_init = u_0
        v_init = v_0
        # the force in the spring
        f_spring_init = STIFFNESS*(u_0-u_p)
        # the unbalanced force, initially equal to the external load increment
        residual = delta_p_eff
        for j in range(MAXITERS):
            # the disp increment according to the initial stiffness
            u_resid = residual/k_eff
            u_next = u_init + u_resid
            v_next = v_init + \
                     G_B*u_resid/STEP - G_B*v_init + STEP*(1.0-0.5*G_B)*a_0
            # we are interested in the total plastic elongation
            oops, u_p = kt_and_up(u_next, v_next, u_p)
            # because we need the spring force at the end of the time step
            f_spring_next = STIFFNESS*(u_next-u_p)
            # so that we can compute the fraction of the incremental force
            # that's equilibrated at the end of the time step
            f_equil = (f_spring_next-f_spring_init) + (k_eff-k_0)*u_resid
            # and finally the incremental forces unbalanced at the end
            # of the time step
            residual = residual-f_equil
            # finish updating the system state
            u_init = u_next
            v_init = v_next
            f_spring_init = f_spring_next
            # if the unbalanced load is small enough (the criteria
            # used in practical (MDOFs) programs are energy based)
            # exit the loop - note that we have no
            # plasticization/unloading residual==0 at the end of the
            # first iteration
            if abs(residual)<F_YIELD*1E-6 :
                break
            else:
                print "#", i, j, u_next, residual/F_YIELD
        # now the load increment is balanced by the spring force and
        # increments in inertial and damping forces, we need to compute
        # the full state at the end of the time step, and to change all
        # denominations to reflect the fact that we are starting a new
        # time step.
        v_0 = v_0 + G_B/STEP*(u_init-u_0)-G_B*v_0+STEP*(1-G_B/2)*a_0
        u_0 = u_init
        k_0, u_p = kt_and_up(u_0, v_0, u_p)
        t_0 = t_1
        p_0 = p_1
        a_0 = ( p_0 - DAMPING*v_0 - STIFFNESS*(u_0-u_p) )/MASS

if __name__ == '__main__':
    main()
