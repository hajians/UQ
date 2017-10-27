#! /usr/bin/env python2.7

class MCMC(object):
    '''
    MCMC implements Monte-Carlo-Markov-Chain algorithms
    '''

    def __init__(self, Density, Proposal, Drawing_func, initial, accept=0.5):
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

        self.accept = accept

        self.num_proposals = 0
        
    def run(self, max_iter = 2000, burning=200):
        '''
        run MCMC.
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
                
            if ( (alpha > self.accept) & (it > burning) ):
                x_old = x_prop
                self.density_samples.append(x_prop)

            if it%100==1:
                print "iter: ", it, ", # samples: ", len(self.density_samples)

            self.num_proposals += 1
            
            
    def write(self, filename):
        '''
        write the samples into the file
        '''

        with open(filename, "w") as output:
            for sample in self.density_samples:
                for coef in sample[:-1]:
                    output.write(str(coef)+",")
                output.write(str(sample[-1])+"\n")
