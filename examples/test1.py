#! /usr/bin/env python2.7

from UQuant.SemilinearSystem import SemiLinSystem

# build a pipe
pipe = SemiLinSystem(1.0, 5.0, 0.0, 1.0, 0.005, 7, 0.05)

vec_c = [0.05, 0.01, 0.0, 0.04, 0.0, 0.005, 0.005]

pipe.run(vec_c)
# get the info of the pipe
pipe.info()

pipe.get_lambda_average(vec_c)
pipe.get_presure_drop(time_instance=13)

