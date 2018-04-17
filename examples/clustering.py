#! /usr/bin/env python2.7

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AffinityPropagation, MeanShift, SpectralClustering

from UQuant.SemilinearSystem import SemiLinSystem

## reading data
# "results/samples-11-v0.0.dat"
df = pd.read_csv("results/uq_pcn.dat", header=None)

threshold = 0.75

df = df[ df.iloc[:,-1] > df.iloc[:,-1].max()*threshold ].iloc[:,:-1]

## physical settings
c_sound = 1.0
t_final = 5.0
x_l, x_r = [0.0, 1.0]
dx = 0.005
boundary_eps = 0.05
expan_coef = len(df.columns)
pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)

## true friction coefficient
true_friction = [0.185938, -0.0519335, 0., 0., -0.0696583, 0.0336323, 0., 0., \
                 0.0348292, -0.0121076, 0., 0., -0.00773981, 0.00105987, 0., 0., 0., \
                 -0.000641154, 0., 0., -0.00278633, 0.00250158, 0., 0., 0.00386991, \
                 -0.00179107, 0., 0., -0.0014216, 0.000230816, 0., 0., 0., \
                 -0.000179701, 0., 0., -0.000859979, 0.000838478, 0., 0., 0.00139317]

pipe_true = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, len(true_friction), boundary_eps)

time_ins = 20

## clustering
X = df.values

### KMeans model
n_clusters = 4
model = KMeans(n_clusters=n_clusters, random_state=100)
#model = SpectralClustering(n_clusters=n_clusters, n_jobs=2, random_state=0)
### AffinityPropagation model
#model = AffinityPropagation()

model.fit(X)

if isinstance(model, KMeans):
    pred = model.predict(X)    
elif isinstance(model, SpectralClustering):
    pred = model.labels_
elif isinstance(model, AffinityPropagation):
    n_clusters = len(model.cluster_centers_indices_)
    pred = model.predict(X)


mean_cluster = []
clusters = []
fig = plt.figure(figsize=(10, 5))

show_clusters = range(n_clusters)
n_show_clusters = len(show_clusters)

for idx, cluster in enumerate(show_clusters):

    ax = fig.add_subplot(1,n_show_clusters,idx+1)
    
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18, rotation=25)
    plt.xlabel("$x$", fontsize=24)
    
    mean_cluster.append( df.iloc[pred==cluster, :].mean().values )
    pipe.get_lambda_average(mean_cluster[idx])
    ax.plot(pipe.mesh, pipe.lambda_avg,
             label="$E(\Lambda_"+str(idx+1)+")(x)$")

    pipe_true.get_lambda_average(true_friction)
    ax.plot(pipe_true.mesh, pipe_true.lambda_avg, "--",
             label="$\lambda_{true}(x)$")

    #ax.legend(loc="upper left", prop={'size': 16}, frameon=False)
    ax.set_aspect(aspect=6)
    ax.legend(loc='upper left', bbox_to_anchor=(0.05, 1.275), fontsize=12)
              

    if cluster!=0:
        ax.tick_params(axis='y', labelleft=False)

plt.tight_layout()
plt.savefig("results/cluster_13.pgf")
plt.show(block=False)

# plot pressure drop

fig = plt.figure(figsize=(12, 6))

pipe_true.run(true_friction)
pipe_true.get_presure_drop(time_instance=time_ins)


markers = ['x', '^', 'P', 'o', '*', 'p']

for cluster in range(n_clusters):
    ax = fig.add_subplot(1,n_clusters,cluster+1)
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    plt.xlabel("$t$", fontsize=24)
    
    mean_cluster.append( df.iloc[pred==cluster, :].mean().values )
    pipe.run(mean_cluster[cluster])
    pipe.get_presure_drop(time_instance=time_ins)
    ax.plot(pipe.timeslices, pipe.pressure_drop, markers[cluster % len(markers)],
             markersize=10)
    ax.set_title("$\Lambda_"+str(cluster+1)+"$",
                 fontsize=24)
    ax.plot(pipe_true.timeslices, pipe_true.pressure_drop,
            "--", label="$\delta p_{true}^n$")

    if cluster!=0:
        ax.tick_params(axis='y', labelleft=False)
    else:
        plt.ylabel("$\delta p_h^n$", fontsize=24)
        
plt.legend(loc=5, borderaxespad=0.0, prop={'size': 20}, frameon=False)
plt.savefig("results/pressuredrop_cluster_13.pgf")
plt.show()

