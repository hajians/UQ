#! /usr/bin/env python2.7

class MCMC(object):
    '''
    MCMC implements Monte-Carlo-Markov-Chain algorithms
    '''

    def __init__(self, Density, Proposal, Drawing_func, initial):
        '''
        Initialize an MCMC algorithm.

        Density : density function whose normalization constant is not known

        Proposal: proposal function

        Drawing_func: a function that draws sample from Proposal
        
        '''
        self.density = Density

        self.proposal = Proposal

        self.sample_proposal = Drawing_func

        self.initial = initial

        self.density_samples = []

        self.accept = 0.5

    def run(self, max_iter = 2000, burning=200):
        '''
        run MCMC.
        '''

        x_old = self.initial
        
        for it in range(max_iter):
            # draw from proposal
            x_prop = self.sample_proposal(x_old)

            # compute alpha !!! assumption: symmetric proposal
            alpha = min(1.0, self.density(x_prop)/self.density(x_old) )

            if ( (alpha > self.accept) & (it > burning) ):
                x_old = x_prop
                self.density_samples.append(x_prop)

        
            if it%100==1:
                print "iter: ", it
