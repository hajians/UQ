#! /usr/bin/env python2.7

from ctypes import cdll, c_double, c_bool, c_char, c_int, POINTER
import numpy as np
import os

_DIR = os.path.dirname(os.path.abspath(__file__))
lib = cdll.LoadLibrary(_DIR+'/lib/CWrapper.so')

class SemiLinSystem(object):

    def __init__(self, c_sound, t_final, x_l, x_r, dx, lambda_len, eps = 10**-2):
        '''
        Initialize a SemiLinearSystem associated to a pipe.

        Parameters
        ----------

        c_sound: float
        sound speed

        t_final: float
        final time

        x_l: float
        left boundary of the domain
        
        x_r: float
        right boundary of the domain

        dx: float
        mesh size

        lambda_len: int
        expansion number for the friction coefficients

        eps: float
        a number for computing pressure drop at both ends of the pipe
        
        '''
        self.obj = lib.CSemiLinSystem(c_double(c_sound),
                                      c_double(t_final),
                                      c_double(x_l), c_double(x_r),
                                      c_double(dx), c_int(lambda_len),
                                      c_double(eps) )

        lib.CNumberofCells.restype = c_int
        self.NumberofCells = lib.CNumberofCells(self.obj)
        
        self.lambda_len = lambda_len
        
        self.boundary_p = None
        self.timeslices = None
        self.pressure_drop = None

        self.mesh = [x_l + i*dx for i in range(self.NumberofCells)]
        
    def info(self):
        '''
        Get the info of the pipe.
        '''
        lib.CInfo(self.obj)

    def run(self, coef, write_bool=False):
        '''
        Compute the solution till the given end time.

        Parameters
        ----------

        write_bool: bool
        if True, then the computed solution at each time step is saved in 'output' file.
        
        '''
        if self.lambda_len==len(coef):
            lib.CRun(self.obj, (c_double * len(coef))(*coef), c_bool(write_bool) )
        else:
            print "Error in run: coef does not have correct size"

    def get_presure_drop(self, time_instance=10, inplace=False):
        '''
        computes pressure drop at certain times.

        Parameters
        ----------
        
        time_instance: int
        pressure drop will be computed at #time_instance uniformly in 
        the interval [0, t_final]

        inplace: bool
        if True, it returns a vector containing pressure drops
        '''

        length = self._CurrentTimeIndex()
        step   = length / time_instance
        idx = range(0,length,step)
        
        self.get_boundary_p()
        
        self.pressure_drop = self.boundary_p[idx,0] - self.boundary_p[idx,1]
        self.timeslices = self._TimeSlices()[idx]

        if inplace==False:
            return self.pressure_drop
        
    def get_boundary_p(self):
        '''
        Returns the boundary values of P into self.boundary_p
        '''
        length = self._CurrentTimeIndex() + 1

        self.boundary_p = np.empty([length,2], dtype=c_double)
        
        lib.CBoundaryValueP_Left.restype = POINTER(c_double)
        lib.CBoundaryValueP_Right.restype = POINTER(c_double)

        data_left = lib.CBoundaryValueP_Left(self.obj)
        data_right = lib.CBoundaryValueP_Right(self.obj)
        
        for i in range(0,length):
            self.boundary_p[i,:] = [data_left[i], data_right[i]]

    def get_lambda_average(self, vec_coef):
        '''
        Gets the lambda average array evaluated on the mesh for the 
        given friction function.

        Parameters
        ----------

        vec_coef: the vector representing the friction coefficient

        Returns
        -------

        self.lambda_avg: numpy.array(double)
        nodal values of the friction function evaluated on the mesh
        
        '''

        self.lambda_avg = np.empty(self.NumberofCells, dtype=c_double)

        lib.CLambda_Average.restype = POINTER(c_double)

        for i in range(self.NumberofCells):
            self.lambda_avg[i] = lib.CLambda_Average(self.obj, (c_double * len(vec_coef))(*vec_coef))[i+1]

    def get_current_lambda_average(self):
        '''
        Get the lambda average array. This is updated after each self.run command.
        '''

        self.lambda_avg = np.empty(self.NumberofCells, dtype=c_double)

        lib.CGetLambda_Average.restype = POINTER(c_double)

        for i in range(self.NumberofCells):
            self.lambda_avg[i] = lib.CGetLambda_Average(self.obj)[i+1]

        
        
    def _TimeSlices(self):
        '''
        Get the time slices, i.e., a vector of type [0, dt, 2*dt, ..., T]
        '''
        length = self._CurrentTimeIndex() + 1
        
        lib.CTimeSlices.restype = POINTER(c_double)

        out = np.empty([length,1], dtype=c_double)
        for i in range(0,length):
            out[i] = lib.CTimeSlices(self.obj)[i]

        return out
    
    def _Write2File(self, filename, append):
        '''
        Write the solution at the current time to a file.
        '''
        lib.CWrite2File(self.obj, filename, c_bool(append))

    def _CurrentTimeIndex(self):
        '''
        Returns the current time index.
        '''
        return lib.CCurrentTimeIndex(self.obj)

    
    
if __name__ == "__main__":

    import matplotlib.pyplot as plt

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
