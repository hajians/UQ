#! /usr/bin/env python2.7

from UQuant.mcmc import MCMC
from UQuant.SemilinearSystem import SemiLinSystem
from numpy.random import normal
from numpy import pi, exp, dot
from numpy import empty
from math  import isnan
from gc    import collect

import pandas as pd

import matplotlib.pyplot as plt


# build a pipe

## physical settings
c_sound = 1.0
t_final = 5.0
x_l, x_r = [0.0, 1.0]
dx = 0.005
boundary_eps = 0.05

## true friction coefficient
true_friction = [0.185938, -0.0519335, 0., 0., -0.0696583, 0.0336323, 0., 0., \
                 0.0348292, -0.0121076, 0., 0., -0.00773981, 0.00105987, 0., 0., 0., \
                 -0.000641154, 0., 0., -0.00278633, 0.00250158, 0., 0., 0.00386991, \
                 -0.00179107, 0., 0., -0.0014216, 0.000230816, 0., 0., 0., \
                 -0.000179701, 0., 0., -0.000859979, 0.000838478, 0., 0., 0.00139317]

true_expan_coef = len(true_friction)

time_ins = 40

## construct and run the true pipe
pipe_true = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, true_expan_coef, boundary_eps)
pipe_true.run(true_friction)
# y_obs = normal(0.0, 0.05, time_ins) + \
#         pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)
y_obs = pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)

# construct a pipe for computation
## stochastic settings
# uni_prior_down = [0.0, -0.025, -0.025, -0.025, -0.025, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01]
# uni_prior_up   = [0.45, 0.025, 0.025, 0.025, 0.025, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]

uni_prior_down = [0.0, -0.025, -0.025, -0.025, -0.025, -0.01, -0.01, -0.01, -0.01]
uni_prior_up   = [0.45, 0.025, 0.025, 0.025, 0.025, 0.01, 0.01, 0.01, 0.01]


sigma_normal   = 0.001

# initial_point_mcmc = [0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# initial_point_mcmc = [0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
initial_point_mcmc = [ 0.18562064, -0.02437987, -0.01683006,  0.0036621 , -0.02162891, \
                       -0.00511077, -0.0002852 , -0.00790147,  0.00633092]

expan_coef = len(initial_point_mcmc)

pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)


def uni_prior(x):
    '''
    uniform prior

    Parameters
    ----------

    x: float * len(x)
    '''

    volume = 1.0

    # for i in range(len(x)):
    #     volume *= abs(uni_prior_up[i] - uni_prior_down[i])
    # *** do not normalize: in high-dim we divide by very small numbers ***

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
    # for i in x:
    #     if i<0: return 0.0
     
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

    filename = "samples-9-v8.3.dat"

    # df = pd.read_csv("tmp.dat", header=None)
    # initial_point_mcmc = df.iloc[-1,0:-1].values
    linestyles = ['-.', '-', '--', ':']
    counter = 0
    plt.style.use('grayscale')

    for i in range(1,len(true_friction)+1,2):
        tmp = true_friction[:i]
        pipe_tmp = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, len(tmp), boundary_eps)

        # get the info of the pipe
        # pipe_tmp.info()
        
        pipe_tmp.get_lambda_average(tmp)
        
        #plt.plot(pipe.timeslices, pipe.pressure_drop, "o-")
        if i in [1,9,11,19,25,41]:
            plt.plot(pipe_tmp.mesh, pipe_tmp.lambda_avg, 
                     linestyle=linestyles[counter%len(linestyles)],
                     label="$N="+str(i/2)+"$")
            counter += 1
    

    plt.legend(loc="best", prop={'size': 16}, frameon=False)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.xlabel("$x$", fontsize=24)
    plt.ylabel("$\lambda^N(x)$", fontsize=24)
    plt.tight_layout()
    plt.savefig("results/figure_fric_bump.pgf")
    plt.show()
    
    # clean memory
    collect()

    pipe.info()
    ### instantiate an MCMC sampler
    mcmc = MCMC(density, proposal_density, draw_from_proposal, initial_point_mcmc)

    # run the MCMC sample
    mcmc.run(max_iter = 200000, burning=5000)
    # write the samples into a file
    mcmc.write(filename, write_prob=True)
