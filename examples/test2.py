#! /usr/bin/python2.7

'''A script to run UQ for a pipe.

Parameters
----------

expan_coef: int
size of friction coefficient vector

uni_prior_down: float * expan_coef
lower bound for the uniform prior

uni_prior_up: float * expan_coef
upper bound for the uniform prior

initial_point_mcmc: float * expan_coef
initial point in the MCMC sampler

c_sound: float
sound speed

t_final: float
final time

x_l: float
left boundary of the domain

x_r: float
right boundary of the domain

dx: float
mesh size

boundary_eps: float
a number for computing pressure drop at both ends of the pipe

pipe_true: SemiLinSystem
an object that represents and contains attributes of the true pipe

pipe: SemiLinSystem 
an object that represents and contains attributes
of a pipe used during MCMC for sampling

'''

from UQuant.mcmc import MCMC
from UQuant.SemilinearSystem import SemiLinSystem
from numpy.random import normal
from numpy import pi, exp, dot
from numpy import empty
from math  import isnan

# stochastic settings
uni_prior_down = [0.0]
uni_prior_up   = [0.5]

sigma_normal   = 0.05
initial_point_mcmc = [0.45]
expan_coef = 1

# physical settings
c_sound = 1.0
t_final = 5.0
x_l, x_r = [0.0, 1.0]
dx = 0.005
boundary_eps = 0.05

# true friction coefficient
true_friction = [0.075]
time_ins = 20

# construct and run the true pipe
pipe_true = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)
pipe_true.run(true_friction)
y_obs = normal(0.0, 0.1, time_ins) + \
        pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)

# construct a pipe for computation
pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)

def uni_prior(x):
    '''
    uniform prior

    Parameters
    ----------

    x: float * len(x)
    '''

    volume = 1.0

    for i in range(len(x)):
        volume *= abs(uni_prior_up[i] - uni_prior_down[i])

    normalize = 1.0 / volume

    for i in range(len(x)):
        if ( ((x[i]<uni_prior_up[i]) & (x[i]>uni_prior_down[i]))==False ):
            return 0.0

    return normalize

def proposal_density(new, old):
    '''
    proposal density function

    Parameters
    ----------
    new: float * len(new)

    old: float * len(old)
    '''

    length = len(new)

    out = empty(length, dtype=float)

    for i in range(length):
        out[i] = 1.0/(2*pi*sigma_normal**2) * \
                 exp( -0.5*(new[i] - old[i])**2 / sigma_normal**2 )
    
    return out

def draw_from_proposal(old):
    '''
    draw from the proposal

    Parameters
    ----------
    old: float * len(old)
    '''
    return normal(old, sigma_normal, size=len(old))

def likelihood(x):
    '''
    computes the likelihood function

    Parameters
    ----------
    x: float * len(x)
    '''
    # check 
    for i in x:
        if i<0: return 0.0
     
    pipe.run(x)
    
    # check if we have negative friction
    pipe.get_current_lambda_average()
    for i in pipe.lambda_avg:
        if i<0: return 0.0
    
    S = pipe.get_presure_drop(time_instance=time_ins, inplace=False)

    out = exp( -0.5 * dot(S-y_obs, S-y_obs) )

    return out
    
def density(x):
    '''
    density function for the MCMC: 
    multiplication of the likelihood and the prior distribution

    Parameters
    ----------
    x: float * len(x)

    '''
    PRIOR = uni_prior(x)
    if PRIOR > 10.0**-8:
        return likelihood(x) * PRIOR
    else:
        return 0.0

if __name__ == "__main__":

    ### instantiate an MCMC sampler
    mcmc = MCMC(density, proposal_density, draw_from_proposal, initial_point_mcmc)

    # run the MCMC sample
    mcmc.run(max_iter = 2000, burning=1000)
    # write the samples into a file
    mcmc.write("samples.dat")


    
