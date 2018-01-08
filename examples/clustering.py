#! /usr/bin/env python2.7

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AffinityPropagation, MeanShift, SpectralClustering

from UQuant.SemilinearSystem import SemiLinSystem

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


## clustering
threshold = 0.95

df = pd.read_csv("samples-9-v7.0.dat", header=None)
df = df[ df.iloc[:,-1] > df.iloc[:,-1].max()*threshold ].iloc[:,:-1]

X = df.values

### KMeans model
n_clusters = 3
model = KMeans(n_clusters=n_clusters, random_state=0)
model = SpectralClustering(n_clusters=n_clusters)
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

for cluster in range(n_clusters):
    mean_cluster.append( df.iloc[pred==cluster, :].mean().values )
    pipe.get_lambda_average(mean_cluster[cluster])
    plt.plot(pipe.mesh, pipe.lambda_avg)

pipe_true.get_lambda_average(true_friction)
plt.plot(pipe_true.mesh, pipe_true.lambda_avg, "--")
plt.show()
