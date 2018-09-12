from UQuant.SemilinearSystem import SemiLinSystem

from numpy.random import normal
from numpy import pi, exp, dot, sqrt
from numpy import empty

from scipy.stats import norm as scipy_norm
from scipy.special import erfinv

from math  import isnan

def isclose(a, b, rel_tol=1e-06, abs_tol=0.0):
    '''
    check if two floats are close
    '''
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def iscloselist(a, List):
    ''' 
    check if a float is close to an element of a list
    '''
    for idx, item in enumerate(List):
        if isclose(a, item):
            return idx
    return None

class BayesianInversion(object):
    '''
    A template class for Bayesian Inversion.
    '''
    def __init__(self, c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps, time_ins):
        from UQuant.SemilinearSystem import SemiLinSystem
    
        self.c_sound, self.t_final = c_sound, t_final
        
        self.x_l, self.x_r = x_l, x_r
        
        self.dx, self.expan_coef, self.boundary_eps = dx, expan_coef, boundary_eps
        
        self.time_ins = time_ins
        
        self.pipe = SemiLinSystem(c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps)

    @staticmethod
    def EvalY_obs(timeslices, y_obs, y_obs_times):
        '''
        evaluates a given time_instances on y_obs.

        Parameters
        ----------

        time_instances: float * len(time_instances)

        '''

        out = empty(len(timeslices), dtype=float)    

        for idx, time in enumerate(timeslices):
            inList = iscloselist(time, y_obs_times)
            if inList!=None:
                out[idx] = y_obs[inList]
            else:
                right = bisect(y_obs_times, time)
                left  = right - 1
                if left < 0:
                    out[idx] = y_obs[right]
                elif right >= len(y_obs_times):
                    out[idx] = y_obs[left]
                else:
                    out[idx] = 0.5*(y_obs[right]+y_obs[left])
        return out
    
class MCMCBayesianInversion(BayesianInversion):
    '''
    GasBayesianInversion solves a Bayesian Inverse 
    problem for Gas pipe friction coefficients.
    '''
    
    def __init__(self, c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps, time_ins):
        BayesianInversion.__init__(self, c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps, time_ins)
         
    
    def uni_prior(self, x):
        '''
        uniform prior

        Parameters
        ----------

        x: float * len(x)
        '''

        uni_prior_up = self.uni_prior_range["up"]
        uni_prior_down = self.uni_prior_range["down"]
        
        volume = 1.0

        for i in range(len(x)):
            volume *= abs(uni_prior_up[i] - uni_prior_down[i])

        normalize = 1.0 / volume

        for i in range(len(x)):
            if ( ((x[i]<uni_prior_up[i]) & (x[i]>uni_prior_down[i]))==False ):
                return 0.0

        return normalize

    def proposal_density(self, new, old):
        '''
        proposal density function

        Parameters
        ----------
        new: float * len(new)

        old: float * len(old)
        '''

        length = len(new)

        out = empty(length, dtype=float)

        for i in range(length):
            out[i] = 1.0/(2*pi*self.sigma_normal**2) * \
                     exp( -0.5*(new[i] - old[i])**2 / self.sigma_normal**2 )

        return out

    def draw_from_proposal(self, old):
        '''
        draw from the proposal

        Parameters
        ----------
        old: float * len(old)
        '''
        return normal(old, self.sigma_normal, size=len(old))

    def likelihood(self, x):
        '''
        computes the likelihood function

        Parameters
        ----------

        pipe: SemiLinSystem

        x: float * len(x)
        '''

        pipe = self.pipe
        
        pipe.run(x, progress_bool=False)

        # check if we have negative friction
        pipe.get_current_lambda_average()
        for i in pipe.lambda_avg:
            if i<0.0: return 0.0

        S = pipe.get_presure_drop(time_instance=self.time_ins, inplace=False)

        Proj_y_obs = self.EvalY_obs(pipe.timeslices, self.y_obs, self.y_obs_times)

        out = exp( -0.5 * dot(S - Proj_y_obs, S - Proj_y_obs) / self.epsilon )

        #print (out, self.epsilon)
        
        return out

    def density(self, x):
        '''
        density function for the MCMC: 
        multiplication of the likelihood and the prior distribution

        Parameters
        ----------
        x: float * len(x)

        '''
        PRIOR = self.uni_prior(x)
        if PRIOR > 10.0**-8:
            return self.likelihood(x) * PRIOR
        else:
            return 0.0

    def run(self, conf, y_obs, y_obs_times, max_iter=2000, burning=500, jupyter=False):
        '''
        run_mcmc takes an epsilon and compute the y_obs and then run mcmc using those observations.

        Returns
        -------

        samples: list
        '''
        from UQuant.mcmc import MCMC
              
        self.epsilon = conf["epsilon"]
        self.sigma_normal = conf["sigma_normal"]
        self.initial_point_mcmc = conf["initial_point_mcmc"]
        self.uni_prior_range = conf["uni_prior_range"]
        
        self.y_obs = y_obs
        self.y_obs_times = y_obs_times
            
        self.sampler = MCMC(self.density, 
                            self.proposal_density, 
                            self.draw_from_proposal, 
                            self.initial_point_mcmc)
        
        
        self.sampler.run(max_iter=max_iter, burning=burning, jupyter=jupyter)
        
class PCNBayesianInversion(BayesianInversion):
    '''
    GasBayesianInversion solves a Bayesian Inverse 
    problem for Gas pipe friction coefficients.
    '''
    
    def __init__(self, c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps, time_ins):
        BayesianInversion.__init__(self, c_sound, t_final, x_l, x_r, dx, expan_coef, boundary_eps, time_ins)       
    
    def normal_to_uniform(self, x):
        """normal_to_uniform converts a given vector x with normal
        distribution (0, I) into a vector with uniform distribution.

        x: [float] * len(x)

        returns: [float] * len(x)
        """

        uni_prior_up = self.uni_prior_range["up"]
        uni_prior_down = self.uni_prior_range["down"]
        
        length = len(x)

        output = scipy_norm.cdf(x)

        for i in range(length):
            limit = (uni_prior_up[i] - uni_prior_down[i])
            output[i] = output[i]*limit + uni_prior_down[i]

        return output

    def uniform_to_normal(self, x):
        """converts a uniform random variable to a normal distribution
        """
        import copy
        
        uni_prior_up = self.uni_prior_range["up"]
        uni_prior_down = self.uni_prior_range["down"]        

        length = len(x)

        output = copy.deepcopy(x)

        for i in range(length):
            limit = (uni_prior_up[i] - uni_prior_down[i])
            tmp   = (output[i] - uni_prior_down[i]) / limit
            output[i] = sqrt(2.0) * erfinv(2.0*tmp - 1.0)
        return output

    def Nlikelihood(self, x):
        '''computes the negative of likelihood function

        Parameters
        ----------
        x: [float] * len(x)

        here x is supposed to have a normal distribution with zero mean
        and standard deviation of size 1.
        '''

        LARGE_NUM = 1000000.0

        pipe = self.pipe
        
        pipe.run(x, progress_bool=False)

        # check if we have negative friction
        pipe.get_current_lambda_average()
        for i in pipe.lambda_avg:
            if i<0.0: return LARGE_NUM
      
        S = pipe.get_presure_drop(time_instance=self.time_ins, inplace=False)
        
        Proj_y_obs = self.EvalY_obs(pipe.timeslices, self.y_obs, self.y_obs_times)            
        
        out = 0.5 * dot(S - Proj_y_obs, S - Proj_y_obs) / self.epsilon

        return out

    def run(self, conf, y_obs, y_obs_times, max_iter=2000, burning=500, jupyter=False):
        '''
        run_mcmc takes an epsilon and compute the y_obs and then run mcmc using those observations.

        Returns
        -------

        samples: list
        '''
        from UQuant.pcn import PCN
              
        self.epsilon = conf["epsilon"]
        self.sigma_normal = conf["sigma_normal"]
        self.initial_point_mcmc = conf["initial_point_mcmc"]
        self.uni_prior_range = conf["uni_prior_range"]
        
        self.y_obs = y_obs
        self.y_obs_times = y_obs_times
            
        self.sampler = PCN(self.Nlikelihood, 
                           self.uniform_to_normal, 
                           self.normal_to_uniform, 
                           self.initial_point_mcmc)

        beta = conf["beta"]
        
        self.sampler.run(max_iter=max_iter, burning=burning, beta=beta, jupyter=jupyter)
