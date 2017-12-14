#! /usr/bin/python2.7

'''A script that contains data concerning the pipe, discretization,
prior- and likelihood functions. In this script the friction
coefficient is assumed to be scalar.


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

from UQuant.SemilinearSystem import SemiLinSystem
from UQuant.mcmc import MCMC
from numpy.random import normal
from numpy import pi, exp, dot
from numpy import empty
from bisect import bisect
from math  import isnan

# stochastic settings
uni_prior_down = [0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0]
uni_prior_up   = [0.5, 0.05,0.05, 0.05, 0.05, 0.05, 0.05]

sigma_normal   = 0.05
initial_point_mcmc = [0.45, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]
expan_coef = 7

# physical settings
c_sound = 1.0
t_final = 5.0
x_l, x_r = [0.0, 1.0]
dx = 0.005
boundary_eps = 0.05

# true friction coefficient
true_friction = [0.075, 0.015, 0.035, 0.025, 0.035, 0.001, 0.03]
time_ins = 20

# construct and run the true pipe
pipe_true = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)
pipe_true.run(true_friction)
y_obs = normal(0.0, 0.1, time_ins) + \
        pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)
y_obs_times = pipe_true.timeslices


def isclose(a, b, rel_tol=1e-06, abs_tol=0.0):
    '''
    check if two floats are close
    '''
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def iscloselist(a, List):
    ''' 
    check if a float is close to an element of a list
    '''
    for idx, item in enumerate(List):
        if isclose(a, item):
            return idx
    return None

def EvalY_obs(timeslices):
    '''
    evaluates a given time_instances on y_obs.

    Parameters
    ----------
    
    time_instances: float * len(time_instances)

    '''

    out = empty(len(timeslices), dtype=float)    

    for idx, time in enumerate(timeslices):
        inList = iscloselist(time,y_obs_times)
        if inList!=None:
            out[idx] = y_obs[inList]
        else:
            right = bisect(y_obs_times, time)
            left  = right - 1
            if left < 0:
                out[idx] = y_obs[right]
            elif right >= len(y_obs_times):
                out[idx] = y_obs[left]
            else:
                out[idx] = 0.5*(y_obs[right]+y_obs[left])
    return out


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

    pipe: SemiLinSystem

    x: float * len(x)
    '''
    # check 
    for i in x:
        if i<0: return 0.0
    
    pipe.run(x, progress_bool=False)
    
    # check if we have negative friction
    pipe.get_current_lambda_average()
    for i in pipe.lambda_avg:
        if i<0: return 0.0
    
    S = pipe.get_presure_drop(time_instance=time_ins, inplace=False)
    
    Proj_y_obs = EvalY_obs(pipe.timeslices)

    out = exp( -0.5 * dot(S - Proj_y_obs, S - Proj_y_obs) )

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

if __name__=="__main__":

    import matplotlib.pyplot as plt
    import gc

    # plt.plot(pipe_true.timeslices, pipe_true.pressure_drop, '--', 
    #          label="true")

    ### instantiate an MCMC sampler
    dx = 0.005

    samples = {}

    # pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)
    # pipe.run([0.075])
    # pipe.get_presure_drop(inplace=True)
    # plt.plot(pipe.timeslices, pipe.pressure_drop, "-x")
    # #plt.plot(pipe_true.timeslices, pipe_true.pressure_drop, "-o")
    # #plt.plot(y_obs_times, y_obs, "-o")
    # plt.plot(pipe.timeslices, EvalY_obs(pipe.timeslices), "-*")
    # plt.show(block=False)
    
    
    for i in range(0,11):
        if i%10==1:
            pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)
            pipe.info()
            mcmc = MCMC(density, proposal_density, draw_from_proposal, initial_point_mcmc)
            mcmc.run(2000,burning=500)
            mcmc.write("samples"+str(i))
            samples[dx] = mcmc.density_samples
        dx += 0.0001

        gc.collect()            # collect un-used objects

    # plt.legend(loc=2, borderaxespad=0.0)
    # plt.show(block=False)

    # Samples = {}
    # for idx in samples:
    #     Samples[idx] = []
    #     for point in samples[idx]:
    #         Samples[idx].append(point[0])

    # for idx in Samples:
    #     plt.hist(Samples[idx], bins=40, normed=True, label=str(idx))

    # plt.legend(loc=2)
