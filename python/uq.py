#! /usr/bin/python2.7

'''
Implements UQ for a pipe
'''

from python.mcmc import MCMC
from python.SemilinearSystem import SemiLinSystem
from numpy.random import normal
from numpy import pi, exp, dot
from numpy import empty
from math  import isnan

# stochastic settings
uni_prior_down = [0.0, 0.0,0.0, 0.0, 0.0]
uni_prior_up   = [0.5, 0.05,0.05, 0.05, 0.05]
sigma_normal   = 0.05
initial_point_mcmc = [0.45, 0.04, 0.04, 0.04, 0.04]
expan_coef = 5

# physical settings
c_sound = 1.0
t_final = 5.0
x_l, x_r = [0.0, 1.0]
dx = 0.005
boundary_eps = 0.05

# true friction coefficient
true_friction = [0.075, 0.015, 0.035, 0.025, 0.035]
time_ins = 20

# construct and run the true pipe
pipe_true = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)
pipe_true.run(true_friction)
# y_obs = normal(0.0, 0.1, time_ins) + \
#         pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)
y_obs = pipe_true.get_presure_drop(time_instance=time_ins, inplace=False)

# construct a pipe for computation
pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)

def uni_prior(x):
    '''
    uniform prior
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
    '''
    return normal(old, sigma_normal)

def likelihood(x):

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

    if isnan(out):
        print "Mistake"
        print x
    return out
    
def density(x):
    '''
    density function for the MCMC
    '''
    PRIOR = uni_prior(x)
    if PRIOR > 10.0**-8:
        return likelihood(x) * PRIOR
    else:
        return 0.0

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    
    pipe_true.info()
    plt.plot(pipe_true.timeslices, pipe_true.pressure_drop)
    plt.plot(pipe_true.timeslices, y_obs)
    plt.show()
    
    mcmc = MCMC(density, proposal_density, draw_from_proposal, initial_point_mcmc)
    # print mcmc.density(0.051)

    mcmc.run(max_iter = 10000, burning=500)
    mcmc.write("samples-fric-0.075.dat")

    # plt.hist(mcmc.density_samples, bins=20)
    # plt.show()

    for idx, coef in enumerate(mcmc.density_samples):
        pipe.get_lambda_average(coef)
        plt.plot(pipe.mesh, pipe.lambda_avg)

    pipe_true.get_lambda_average(true_friction)
    plt.plot(pipe_true.mesh, pipe_true.lambda_avg, '--',
             color="black", linewidth=5.0, label="True friction function")

    plt.legend(bbox_to_anchor=(1.0, 1), loc=1, borderaxespad=0.0)
    plt.xlabel("$x$", fontsize=22)
    plt.ylabel("$\lambda(x)$", fontsize=22)
    plt.show(block=False)


# pipes = []
# for i in range(1000):
#     pipes.append(SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)) 
# for pipe in pipes:
#     pipe.run([true_friction])
#     pipe.get_presure_drop(time_instance=time_ins)
#     plt.plot(pipe.timeslices, pipe.pressure_drop, "-o")    

# plt.show(block=False)
    
