#! /usr/bin/env python2.7

from numpy import random
import time

# change the seed
random.seed(int(time.time()))

class MCMC(object):
    '''
    MCMC implements Markov-Chain-Monte-Carlo algorithms.
    '''

    def __init__(self, Density, Proposal, Drawing_func, initial):
        '''
        Initialize an MCMC algorithm.

        Parameters
        ----------
        Density : density function whose normalization constant is not known

        Proposal: proposal function

        Drawing_func: a function that draws sample from Proposal

        initial: the initial sample of the MCMC algorithm
        
        '''
        self.density = Density

        self.proposal = Proposal

        self.sample_proposal = Drawing_func

        self.initial = initial

        self.density_samples = []

        self.num_proposals = 0


    def run(self, max_iter = 2000, burning=200):
        '''run MCMC.

        Parameters
        ----------

        max_iter: int
        
        maximum number of iterations on this call
        
        burning: int 

        ignores the initial samples upto the given number in the
        Markov Chain.

        '''

        if len(self.density_samples)>0:
            x_old = self.density_samples[-1]
            burning = 0
        else:
            x_old = self.initial
        
        for it in range(max_iter):
            # draw from proposal
            x_prop = self.sample_proposal(x_old)

            # compute alpha !!! assumption: symmetric proposal
            try:
                alpha = min(1.0, self.density(x_prop)/self.density(x_old) )
            except:
                print "error in MCMC.run:"
                print self.density(x_prop), x_prop
                print self.density(x_old), x_old
                alpha = 0.0
                
            if ( (alpha > random.uniform(0.0,1.0)) & (it > burning) ):
                x_old = x_prop
                self.density_samples.append(x_prop)

            if it%100==1:
                print "iter: ", it, ", # samples: ", len(self.density_samples)

            self.num_proposals += 1
            
            
    def write(self, filename):
        '''
        write the samples into the file.

        Parameters
        -----------

        filename: str
        
        '''

        with open(filename, "w") as output:
            for sample in self.density_samples:
                for coef in sample[:-1]:
                    output.write(str(coef)+",")
                output.write(str(sample[-1])+"\n")

    
