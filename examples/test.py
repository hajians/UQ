#! /usr/bin/env python2.7

import matplotlib.pyplot as plt

from UQuant.SemilinearSystem import SemiLinSystem

# build a pipe
pipe = SemiLinSystem(1.0, 5.0, 0.0, 1.0, 0.005, 7, 0.05)
# get the info of the pipe
pipe.info()

# check how pressure drop behaves when the friction coefficient is
# updated. Theoretically we have continuous dependence which can
# be seen also numerically.
coef = 0.0
for i in range(0,100):
    vec_c = [0.05 + coef, 0.01, 0.0, 0.04, 0.0, 0.005, 0.005]
    pipe.run(vec_c)
    pipe.get_lambda_average(vec_c)
    pipe.get_presure_drop(time_instance=13)
    if i%10==1:
        plt.plot(pipe.timeslices, pipe.pressure_drop, "-o",
                 label="cf = "+ str(0.05 + coef))
    coef += 0.001

print ">> Computation Done"

plt.legend(loc=2, borderaxespad=0.0)
plt.xticks([round(tn,2) for tn in pipe.timeslices])
plt.xlabel("$t_n$", fontsize=24)
plt.ylabel("$\delta p_h(t_n)$", fontsize=24)
plt.tight_layout()
plt.show(block=True)
