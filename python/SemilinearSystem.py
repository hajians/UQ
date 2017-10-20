#! /usr/bin/env python2.7

from ctypes import cdll, c_double, c_bool, c_char, c_int, POINTER
import numpy as np

lib = cdll.LoadLibrary('lib/CWrapper.so')

class SemiLinSystem(object):

    def __init__(self, c_sound, t_final, x_l, x_r, dx, lambda_len, eps = 10**-2):
        '''
        Initialize a SemiLinearSystem associated to a pipe.
        '''
        self.obj = lib.CSemiLinSystem(c_double(c_sound),
                                      c_double(t_final),
                                      c_double(x_l), c_double(x_r),
                                      c_double(dx), c_int(lambda_len),
                                      c_double(eps) )

        self.boundary_p = None
        self.timeslices = None
        self.pressure_drop = None

    def info(self):
        '''
        Get the info of the pipe.
        '''
        lib.CInfo(self.obj)

    def run(self, coef, write_bool=False):
        '''
        Compute the solution till the given end time.
        '''
        lib.CRun(self.obj, (c_double * len(coef))(*coef), c_bool(write_bool) )

    def get_presure_drop(self, time_instance=10, inplace=False):
        '''
        compute pressure drop at certain times
        '''

        length = self._CurrentTimeIndex() + 1
        step   = length / time_instance
        idx = range(0,length,step)
        
        self.get_boundary_p()
        
        self.pressure_drop = self.boundary_p[idx,0] - self.boundary_p[idx,1]
        self.timeslices = self._TimeSlices()[idx]

        if inplace==False:
            return self.pressure_drop
        
    def get_boundary_p(self):
        '''
        Returns the boundary values of P
        '''
        length = self._CurrentTimeIndex() + 1

        self.boundary_p = np.empty([length,2], dtype=c_double)
        
        lib.CBoundaryValueP_Left.restype = POINTER(c_double)
        lib.CBoundaryValueP_Right.restype = POINTER(c_double)

        data_left = lib.CBoundaryValueP_Left(self.obj)
        data_right = lib.CBoundaryValueP_Right(self.obj)
        
        for i in range(0,length):
            self.boundary_p[i,:] = [data_left[i], data_right[i]]

    def _TimeSlices(self):
        '''
        Get time slices.
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
        Returns current time index
        '''
        return lib.CCurrentTimeIndex(self.obj)

    
    
if __name__ == "__main__":

    import matplotlib.pyplot as plt
    
    pipe = SemiLinSystem(1.0, 5.0, 0.0, 1.0, 0.005, 0, 0.05)
    pipe.info()

    coef = 0.0
    for i in range(500):
        pipe.run([coef])
        pipe.get_presure_drop()
        if i%9==1:
            plt.plot(pipe.timeslices, pipe.pressure_drop,
                     label="cf = "+ str(coef))
        coef += 0.001

    print "computation done"

    plt.legend(loc=2)
    plt.show(block=True)

    pipe.run([coef], True)
    
