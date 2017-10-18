#! /usr/bin/env python2.7

from ctypes import cdll, c_double, c_bool, c_char, c_int, POINTER
import numpy as np

lib = cdll.LoadLibrary('lib/CWrapper.so')

class SemiLinSystem(object):

    def __init__(self, c_sound, t_final, x_l, x_r, dx, lambda_len):
        '''
        Initialize a SemiLinearSystem associated to a pipe.
        '''
        self.obj = lib.CSemiLinSystem(c_double(c_sound),
                                      c_double(t_final),
                                      c_double(x_l), c_double(x_r),
                                      c_double(dx), c_int(lambda_len)
                                      )

        self.boundary_p = None
        
    def info(self):
        '''
        Get the info of the pipe.
        '''
        lib.CInfo(self.obj)

    def run(self, coef):
        '''
        Compute the solution till the given end time.
        '''
        lib.CRun(self.obj, (c_double * len(coef))(*coef) )

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
    
    pipe = SemiLinSystem(1.0, 10.0, 0.0, 1.0, .01, 0)
    pipe.info()
    pipe.run([0.5])

    pipe.get_boundary_p()

    plt.plot(pipe.boundary_p[:,0])
    plt.plot(pipe.boundary_p[:,1])

    plt.show()
