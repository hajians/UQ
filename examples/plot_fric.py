#! /usr/bin/env python2.7

from UQuant.mcmc import MCMC
from UQuant.SemilinearSystem import SemiLinSystem
from numpy.random import normal
from numpy import pi, exp, dot
from numpy import empty
from math  import isnan
from gc    import collect

import matplotlib.pyplot as plt


# build a pipe

## physical settings
c_sound = 1.0
t_final = 5.0
x_l, x_r = [0.0, 1.0]
dx = 0.005
boundary_eps = 0.05

## true friction coefficient
true_friction = [1.01563, 0.0148381, 0., 0., -0.0126651, -0.00960923, 0., 0., \
                 0.00633257, 0.00345932, 0., 0., -0.00140724, -0.000302819, 0., 0., \
                 0., 0.000183187, 0., 0., -0.000506606, -0.000714736, 0., 0., \
                 0.000703619, 0.000511734, 0., 0., -0.000258472, -0.0000659473, 0., \
                 0., 0., 0.0000513431, 0., 0., -0.00015636, -0.000239565, 0., 0., \
                 0.000253303]

true_expan_coef = len(true_friction)

time_ins = 20

## construct and run the true pipe
pipe_true = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, true_expan_coef, boundary_eps)
pipe_true.run(true_friction)
y_obs = normal(0.0, 0.1, time_ins) + \
        pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)



# construct a pipe for computation
## stochastic settings
uni_prior_down = [0.0, -0.5, -0.5, -0.5, -0.5]
uni_prior_up   = [2.0, 0.5, 0.5, 0.5, 0.5]

sigma_normal   = 0.05

initial_point_mcmc = [0.45, 0.0, 0.0, 0.0, 0.0]
expan_coef = len(initial_point_mcmc)

pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)
pipe.info()

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
     
    pipe.run(x, progress_bool=False)
    
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

    
    for i in range(1,len(true_friction),2):
        tmp = true_friction[:i]
        pipe_tmp = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, len(tmp), boundary_eps)

        # get the info of the pipe
        # pipe_tmp.info()
        
        pipe_tmp.get_lambda_average(tmp)
        
        #plt.plot(pipe.timeslices, pipe.pressure_drop, "o-")
        plt.plot(pipe_tmp.mesh, pipe_tmp.lambda_avg)
    
    #plt.show()
    
    # clean memory
    collect()

    
    ### instantiate an MCMC sampler
    mcmc = MCMC(density, proposal_density, draw_from_proposal, initial_point_mcmc)

    # run the MCMC sample
    mcmc.run(max_iter = 2000, burning=200)
    # write the samples into a file
    mcmc.write("samples.dat")


    
