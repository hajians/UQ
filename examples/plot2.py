#! /usr/bin/env python2.7

'''
plot2.py animates the sequence of samples in the Markov chain.
'''

import sys

from numpy.random import normal
from numpy import pi, exp, dot
from numpy import empty

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import animation

from UQuant.SemilinearSystem import SemiLinSystem

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

# construct a pipe for computation
pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)

filename = "results/samples-N3-fric0.075-wNoise.dat"
df = pd.read_csv(filename, header=None)
density_samples = df.values

pipe_true.info()

fig = plt.figure()
ax = plt.axes(xlim=(min(pipe_true.mesh), max(pipe_true.mesh)),
              ylim=(0.0,0.5))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel("$x$", fontsize=24)
plt.ylabel("$\lambda(x)$", fontsize=24)
plt.tight_layout()
text_label = ax.text(1,1, '', transform=ax.transAxes)

lines = [ax.plot([], [], "-", lw=1)[0] for _ in range(len(density_samples)+1)]

def init():
    for line in lines:
        line.set_data([], [])
    return lines

def update(i):
    if i>0:
        pipe.get_lambda_average(density_samples[i])
        for data in pipe.lambda_avg:
            if abs(data)>10**2:
                print data
        lines[i-1].set_data( pipe.mesh, pipe.lambda_avg )
        lines[i-1].set_linewidth(1)
        lines[i-1].set_linestyle("-")
        text_label.set_text('sample: '+str(i))

    pipe_true.get_lambda_average(true_friction)
    lines[i].set_data( pipe_true.mesh, pipe_true.lambda_avg )
    lines[i].set_linewidth(4)
    lines[i].set_linestyle("--")
    
    lines[i-1].set_color( lines[i].get_color() )
    lines[i].set_color("black")
    
    text_label.set_text('sample: '+str(i))

    return lines


if __name__=="__main__":
    anim = animation.FuncAnimation(fig, update, init_func=init,
                                   frames=len(density_samples),
                                   repeat=False, repeat_delay=5000,
                                   blit=True)

    try:
        if '-save'==sys.argv[1]:
            anim.save('UQsamples.gif', dpi=80, writer='imagemagick')
    except:
        print "did not save the GIF file"
        
    plt.show(block=False)


    # plot the mean function, true and the initial sample
    plt.figure()

    pipe.get_lambda_average(df.mean().values)
    Mean = pipe.lambda_avg
    plt.plot(pipe.mesh, Mean, label="$\mathrm{E}(\lambda)(x)$", linewidth=2)
    
    pipe.get_lambda_average(initial_point_mcmc)
    plt.plot(pipe.mesh, pipe.lambda_avg,
             label="$\lambda_{initial}(x)$", linestyle="-.", linewidth=2)
    
    pipe_true.get_lambda_average(true_friction)
    plt.plot(pipe_true.mesh, pipe_true.lambda_avg,
         label="$\lambda_{true}(x)$", linestyle="--", linewidth=2)
    plt.legend(loc=5, borderaxespad=0.0, prop={'size': 20}, frameon=False)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.xlabel("$x$", fontsize=24)
    plt.tight_layout()
    #plt.savefig("../../fig_func_mean.pgf")
    plt.show()

