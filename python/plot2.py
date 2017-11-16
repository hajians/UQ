#! /usr/bin/env python2.7

'''
plot2.py animates the sequence of samples in the Markov chain.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib as mpl
from python.uq import *

filename = "samples-N3-fric0.075-wNoise.dat"
df = pd.read_csv(filename, header=None)

pipe_true.info()

mcmc = MCMC(density, proposal_density, draw_from_proposal, initial_point_mcmc)
mcmc.density_samples = df.values

fig = plt.figure()
ax = plt.axes(xlim=(min(pipe_true.mesh), max(pipe_true.mesh)),
              ylim=(0.0,0.5))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel("$x$", fontsize=24)
plt.ylabel("$\lambda(x)$", fontsize=24)
plt.tight_layout()
text_label = ax.text(1,1, '', transform=ax.transAxes)

lines = [ax.plot([], [], "-", lw=1)[0] for _ in range(len(mcmc.density_samples)+1)]

def init():
    for line in lines:
        line.set_data([], [])
    return lines

def update(i):
    if i>0:
        pipe.get_lambda_average(mcmc.density_samples[i])
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
                                   frames=len(mcmc.density_samples),
                                   repeat=False, repeat_delay=5000,
                                   blit=True)

    try:
        anim.save('UQsamples.gif', dpi=80, writer='imagemagick')
    except:
        print "could not same the GIF file"
    
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

