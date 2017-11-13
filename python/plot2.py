#! /usr/bin/env python2.7

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation

from python.uq import *

filename = "samples-fric-0.075-expan-7.dat"
df = pd.read_csv(filename, header=None)

pipe_true.info()

mcmc = MCMC(density, proposal_density, draw_from_proposal, initial_point_mcmc)
mcmc.density_samples = df.values

fig = plt.figure()
ax = plt.axes(xlim=(min(pipe_true.mesh), max(pipe_true.mesh)),
              ylim=(0.0,0.5))
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

anim = animation.FuncAnimation(fig, update, init_func=init,
                               frames=len(mcmc.density_samples),
                               repeat=True, repeat_delay=5000,
                               blit=True)
plt.show()
