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
        '''
        if self.lambda_len==len(coef):
            lib.CRun(self.obj, (c_double * len(coef))(*coef), c_bool(write_bool) )
        else:
            print "Error in run: coef does not have correct size"

    def get_presure_drop(self, time_instance=10, inplace=False):
        '''
        compute pressure drop at certain times.
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
        Returns the boundary values of P.
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
        Get the lambda average array
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
        Get the time slices.
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

    pipe = SemiLinSystem(1.0, 5.0, 0.0, 1.0, 0.005, 7, 0.05)
    pipe.info()
    # pipe.run([0.05, 0.01, 0.001], write_bool=True)
    coef = [0.05, 0.01, 0.0, 0.04, 0.0, 0.0, 0.0]
    pipe.get_lambda_average(coef)
    
    # plt.plot(range(pipe.NumberofCells), pipe.lambda_avg, "-o")
    # plt.show(block=False)
    
    coef = 0.0
    for i in range(0,300):
        vec_c = [0.05, 0.01, 0.0, 0.04, 0.0, coef, 0.0]
        pipe.run(vec_c)
        pipe.get_lambda_average(vec_c)
        pipe.get_presure_drop()
        if i%10==1:
            plt.plot(pipe.timeslices, pipe.pressure_drop, "-o",
                     label="cf = "+ str(coef))
            # plt.plot(range(pipe.NumberofCells), pipe.lambda_avg, "-")
        coef += 0.001

    # print "computation done"

    plt.legend(bbox_to_anchor=(1.0, 1), loc=2, borderaxespad=0.0)
    plt.show(block=True)

#    pipe.run([coef], True)
