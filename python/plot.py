#! /usr/bin/env python2.7

'''This script plot and animates the time-dependent solution obtained
from the forward solver.

'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation

_jump = 2
_dir  = "post"
_filename_Q = "output_Q"
_filename_P = "output_P"
_ext = ".dat"
_fig_ext = ".png"

_wave_speed = 1

dfQ = pd.read_csv(_filename_Q+_ext, names=["x", "value", "time"])
dfP = pd.read_csv(_filename_P+_ext, names=["x", "value", "time"])

df = dfP
#df["value"] = 0.5 * dfP["value"] - 0.5/_wave_speed * dfQ["value"]

time_slices = df.time.unique()[0::_jump]

fig = plt.figure()
ax = plt.axes(xlim=(df["x"].min(), df["x"].max()),
              ylim=(df["value"].min(), df["value"].max()))
time_text = ax.text(0.05,0.9, '', transform=ax.transAxes)
line, = ax.plot([], [], "-o", lw=1)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def update(i):
    df_loc = df[ df["time"] == time_slices[i] ][["x", "value"]]
    line.set_data( df_loc["x"].values, df_loc["value"].values )
    time_text.set_text('time: '+str(time_slices[i]))

    return line, time_text

anim = animation.FuncAnimation(fig, update, init_func=init,
                               frames=len(time_slices), blit=True)
plt.show()
