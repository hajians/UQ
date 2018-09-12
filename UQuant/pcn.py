#! /usr/bin/python2.7

from numpy import random, exp, sqrt
from numpy.random import normal, uniform
from numpy import empty

import time

# change the seed
random.seed(int(time.time()))

class PCN(object):
    """
    PCN implements a preconditioned Crank-Nicholson algorithm.
    """

    def __init__(self, NLikelihood, dist_to_normal, normal_to_dist, initial):
        '''
        Initialize an MCMC algorithm.

        Parameters
        ----------
        NLikelihood: negative likelihood function

        initial: the initial sample of the MCMC algorithm
        
        '''
        self.Nlikelihood = NLikelihood

        self.dist_to_normal = dist_to_normal

        self.normal_to_dist = normal_to_dist
        
        self.initial = initial

        self.density_samples = []

        self.density_probability = []
        
        self.num_proposals = 0
        
        # stats are used to check how many 
        # proposals are accepted and rejected.
        self.stats = {
            "accepted": 0,
            "rejected": 0,
        }

    @staticmethod
    def sample_proposal(x_old, beta):
        """
        sample_proposal generates a new proposal for the PCN.

        x_old: [float] * len(x_old)
        previous sample.

        beta: float
        constant in PCN proposal step.

        returns:
        x_new: [float] * len(x_old)
        proposal.
        """

        length = len(x_old)
        x_prop  = empty(length, dtype=float)

        for i in range(length):
            x_prop[i] = sqrt(1.0 - beta**2) * x_old[i] + beta * normal(0.0, 1.0)

        return x_prop

    def run(self, max_iter = 2000, burning=200, beta=0.5, jupyter=False):
        """
        run PCN.

        max_iter: int    
        maximum number of iterations on this call
        
        burning: int 
        ignores the initial samples upto the given number in the
        Markov Chain.

        beta: float
        constant in the proposal step of PCN.

        """

        if jupyter:
            from tqdm import tqdm_notebook as tqdm
        else:
            from tqdm import tqdm

        if len(self.density_samples)>0:
            x_old = self.density_samples[-1]
            burning = 0
        else:
            x_old = self.initial

        for it in tqdm(range(max_iter)):
            # draw from proposal
            x_old_normal = self.dist_to_normal(x_old)
            x_prop_normal = self.sample_proposal(x_old_normal, beta)
            x_prop = self.normal_to_dist(x_prop_normal)
            
            # compute alpha 
            try:
                likelihood_old = self.Nlikelihood(x_old)
                likelihood_prop = self.Nlikelihood(x_prop)
                alpha = min(1.0, exp(likelihood_old - likelihood_prop) )
            except:
                print "error in PCN.run:"
                print self.Nlikelihood(x_prop), x_prop
                print self.Nlikelihood(x_old), x_old
                alpha = 0.0
            
            if ( alpha > uniform(0.0,1.0) ):
                x_old = x_prop
                if (it > burning):
                    self.density_samples.append(x_prop)
                    self.density_probability.append(exp(-likelihood_prop))
                    
                    self.stats["accepted"] += 1
            else:
                if (it > burning):
                    self.stats["rejected"] += 1                

            self.num_proposals += 1
                
    def write(self, filename, write_prob=False):
        '''
        write the samples into the file.

        Parameters
        -----------

        filename: str
        
        '''

        with open(filename, "w") as output:
            for i, sample in enumerate(self.density_samples):
                if write_prob:
                    for coef in sample:
                        output.write(str(coef)+",")
                    output.write(str(self.density_probability[i])+"\n")
                else:
                    for coef in sample[:-1]:
                        output.write(str(coef)+",")
                    output.write(str(sample[-1])+"\n")

