#! /usr/bin/python2.7

'''
Implements UQ for a pipe
'''

from python.mcmc import MCMC
from python.SemilinearSystem import SemiLinSystem
from numpy.random import normal
from numpy import pi, exp, dot

# stochastic settings
uni_prior_down = 0.0
uni_prior_up   = 0.5
sigma_normal   = 0.05

# physical settings
c_sound = 1.0
t_final = 5.0
x_l, x_r = [0.0, 1.0]
dx = 0.005
expan_coef = 0
boundary_eps = 0.05

# true friction coefficient
true_friction = 0.05
time_ins = 20

# construct and run the Truth $
pipe_true = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)
pipe_true.run([true_friction])
y_obs = pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)

# construct a pipe for computation
pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)

def uni_prior(x):
    '''
    uniform prior
    '''
    normalize = 1.0 / (uni_prior_up - uni_prior_down)
    
    if ( (x<uni_prior_up) & (x>uni_prior_down) ):
        return normalize
    else:
        return 0.0

def proposal_density(new, old):
    '''
    proposal density function
    '''
    return 1.0/(2*pi*sigma_normal**2) * exp( -0.5*(new-old)**2 / sigma_normal**2 )

def draw_from_proposal(old):
    '''
    draw from the proposal
    '''
    return normal(old, sigma_normal)

def likelihood(x):

    if (x<0):
        return 0.0
    else:
        pipe.run([x])
        S = pipe.get_presure_drop(time_instance=time_ins, inplace=False)

        return exp( -0.5 * dot(S-y_obs, S-y_obs) )
    
def density(x):
    '''
    density function for the MCMC
    '''
    return likelihood(x) * uni_prior(x)    
    

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    
    pipe_true.info()
    
    plt.plot(pipe_true.timeslices, y_obs, "-o")
    plt.show()
    
    mcmc = MCMC(density, proposal_density, draw_from_proposal, 0.06)
    print mcmc.density(0.051)

    mcmc.run(max_iter = 4000, burning=200)

    
