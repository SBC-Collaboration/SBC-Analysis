import numpy as np
from scipy.optimize import minimize as fmin
'''
This module defince the HistFit and GHFfitfun classes, as well as the
GHF_xxxx subclasses of GHFfitfun.

Adding new subclasses to GHFfitfun:
    You can create subclasses to GHFfitfun as needed.  Each new subclass
    should look like the following:
        class GHF<myfunctionname>(GHFfitfun):
            numparams = <number of parameters for this function>
            def calculate(self):
                (calculate self.v, self.dvda, and self.d2vda2, based on
                 self.params and self.binedges)
    Optionally, you can also define pdf(self, xd) and/or cdf(self, xd) -- if
    you don't define them, they'll inherit the methods from GHFfitfun that use
    self.calculate to get the cdf and estimate the pdf
'''

class HistFit:
    ''' A HistFit instance is initialized with the output of a call to
       np.histogram (counts and binedges).  Once created, distributions and
       priors can be added to the fit definition.  The fit is executed by
       calling the .fit() instance method.
       '''
    
    def __init__(self, counts, binedges):
        ''' Initialize an instance with a histogram.  Can be called as:
            hist_out = np.histogram(data)
            fit_obj = HistFit(*hist_out)
            '''
        # the following define the fit to be performed
        self.fitfun = [] # list of distributions included in the fit
        self.paramposts = [0] # list of parameter positions for each dist
        self.params = np.zeros((0), dtype=np.float64) # current parameter set
        self.param_rescale = \
            np.zeros((0), dtype=np.float64) # parameter rescaling for optimizer
        self.prior = [] # list of priors applied to fit parameters
        self.binedges = binedges.copy() # binedges of the histogram to be fit
        self.counts = counts.copy() # counts (integers) in the histogram
        self.cut = self.counts > 0 # flag for bins with zero counts
        self.nd = counts.size # degrees of freedom in the fit
        self.method = 'trust-exact' # optimization method
        
        # the following are outputs created by the fit routine
        self.bestfit = None # best fit parameters
        self.bestfit_err = None # 1-sigma uncertainty on best fit parameters
        self.eigen_vects = None # eigen vectors of the hessian matrix
        self.eigen_vars = None # variances associated with eigen_vects
        self.chisq = None # total chisq
        self.chisq_nd = None # reduced chisq
        self.grad = None # derivatives of chisq w.r.t. parameters
        self.hessian = None # derivatives of chisq w.r.t. parameters
        
        # the following are intended for internal use only
        self.lastfit = None # last set of parameters evaluated
        self.v = None # expected number of counts in each bin
        self.dvda = None # derivatives of v w.r.t. parameters
        self.d2vda2 = None # 2nd derivatves of v w.r.t. parameters
        self.llp = None # log-likelihood contribution from prior
        self.dllpda = None # derivatives of llp w.r.t. parameters
        self.d2llpda = None # 2nd derivatives of llp w.r.t. parameters
        self.chisq_vec = None # Poisson residuals by bin
        self.grad_scaled = None # gradient scaled for optimization
        self.hessian_scale = None # hessian scaled for optimization
        self.fit_out = None # raw output of optimization
    
    def addfun(self, fitfunction, initial_params, rescale = None):
        ''' Add distributions ot the instance before executing the fit.
            fitfunction should be an instance of one of the subclasses of
            GHFfitfun (i.e. issubclass(type(fitfunction),GHFfitfun)) should
            evaluate to True).  Must provide the initial guess for the
            parameters as well, in initial_params.  Optional rescale argument
            rescales parameters for optimizer -- by default, parameters with
            non-zero initial guesses are rescaled by that guess.
            '''
        self.fitfun.append(fitfunction)
        self.paramposts.append(self.paramposts[-1] + fitfunction.numparams_eff)
        self.nd -= fitfunction.numparams_eff
        self.params = np.append(self.params, initial_params)
        fitfunction.update_binedges(self.binedges)
        if rescale is None:
            if type(initial_params) is np.ndarray:
                this_rescale = initial_params.copy()
            else:
                this_rescale = np.float64(initial_params).reshape(-1)
            this_rescale[this_rescale==0] = 1
            self.param_rescale = np.append(self.param_rescale, this_rescale)
        else:
            self.param_rescale = np.append(self.param_rescale, rescale)
    
    def addprior(self, priorfunction):
        ''' Add priors.  Each prior should be callable, should take the full
            parameter list as its sole argument, and should return
            (p, dpda, d2pda2), all numpy arrays with dtype float64.
            p has shape () and is the prior probability of the given parameter
            set.  dpda has shape (n,) for n parameters, and is the first
            derivative of p w.r.t. the parameters.  d2pda2 has shape (n, n) and
            contains the 2nd derivatives of p w.r.t. the parameters.
            '''
        self.prior.append(priorfunction)
    
    def allocate_arrays(self):
        ''' Meant for internal use (called by fit()).  Must be called prior to
            calling GHF_LL(), as it allocates the various numpy arrays used
            in calculating the log-likelihood.
            '''
        self.v = np.zeros((self.counts.size), dtype=np.float64)
        self.dvda = np.zeros((self.params.size, self.counts.size),
                             dtype=np.float64)
        self.d2vda2 = np.zeros((self.params.size, self.params.size,
                                self.counts.size), dtype=np.float64)
        self.dllpda = np.zeros((self.params.size), dtype=np.float64)
        self.d2llpda2 = np.zeros((self.params.size, self.params.size), 
                                dtype=np.float64)
        self.chisq_vec = np.zeros((self.counts.size), dtype=np.float64)
    
    def fit(self):
        ''' Executes the fit, using self.params as the initial guess.  After
            the fit, self.params is the new bestfit, so fit() can be called
            repeatedly if e.g. the optimizer indicates it needed more
            iterations.
            Also calculates uncertainties on fit parameters based on the
            hessian matrix (fisher information matrix) at the best fit point.
            '''
        initial_guess = self.params / self.param_rescale
        self.lastfit = initial_guess + np.nan
        self.allocate_arrays()
        self.fit_out = fmin(self.chisqfun, initial_guess,
                            jac=self.gradfun, hess=self.hessfun,
                            method=self.method)
        self.bestfit = self.fit_out.x * self.param_rescale
        self.params[:] = self.bestfit
        self.GHF_LL()
        (D, self.eigen_vects) = np.linalg.eig(self.hessian)
        self.eigen_vars = 2.0 / D
        self.bestfit_err = np.sqrt((self.eigen_vects**2) @ self.eigen_vars)
    
    def chisqfun(self, newparams):
        ''' Passed to optimizer.  Returns reduced chi-square
            '''
        if not np.all(self.lastfit == newparams):
            np.multiply(newparams, self.param_rescale, out=self.params)
            self.lastfit[:] = newparams
            self.GHF_LL()
        return self.chisq_nd
    
    def gradfun(self, newparams):
        ''' Passed to optimizer.  Returns gradient of reduced chi-square
            w.r.t. scaled parameters
            '''
        if not np.all(self.lastfit == newparams):
            np.multiply(newparams, self.param_rescale, out=self.params)
            self.lastfit[:] = newparams
            self.GHF_LL()
        return self.grad_scaled
    
    def hessfun(self, newparams):
        ''' Passed to optimizer.  Returns hessian of reduced chi-square
            w.r.t. scaled parameters
            '''
        if not np.all(self.lastfit == newparams):
            np.multiply(newparams, self.param_rescale, out=self.params)
            self.lastfit[:] = newparams
            self.GHF_LL()
        return self.hessian_scaled
    
    def GHF_LL(self):
        ''' Log-likelihood calculator.  Adds Poisson log-likelihoods for
            histogram bins with logs of probabilities from each prior.  Updates
            chisq, chisq_nd, grad, and hessian attributes.
            '''
        self.v[:]=0
        self.dvda[:]=0
        self.d2vda2[:]=0
        for i_f in range(len(self.fitfun)):
            self.fitfun[i_f].update_params(self.params[self.paramposts[i_f]:
                                                       self.paramposts[i_f+1]])
            self.v += \
                self.fitfun[i_f].v
            self.dvda[self.paramposts[i_f]:self.paramposts[i_f+1], :] += \
                self.fitfun[i_f].dvda
            self.d2vda2[self.paramposts[i_f]:self.paramposts[i_f+1],
                        self.paramposts[i_f]:self.paramposts[i_f+1], :] += \
                self.fitfun[i_f].d2vda2
        
        self.llp = np.float64(0)
        self.dllpda[:] = 0
        self.d2llpda2[:] = 0
        
        for i_p in range(len(self.prior)):
            (thisp, thisdpda, thisd2pda2) = self.prior[i_p](self.params)
            self.llp += np.log(thisp)
            self.dllpda += thisdpda / thisp
            self.d2llpda2+= (thisd2pda2 / thisp) - \
                (thisdpda * thisdpda[:, np.newaxis] / (thisp**2))
        
        self.chisq_vec[:] = 0
        
        self.chisq_vec[~self.cut] = 2*self.v[~self.cut]
        self.chisq_vec[self.cut] = 2*(self.counts[self.cut] *
                                      np.log(self.counts[self.cut] / 
                                             self.v[self.cut]) +
                                      self.v[self.cut] -
                                      self.counts[self.cut])
        grad_vec = 2 * (1-(self.counts/self.v)) * self.dvda
        
        hess_vec = 2 * (1 - (self.counts/self.v)) * self.d2vda2 + \
            (self.counts/(self.v**2)) * self.dvda * self.dvda[:, np.newaxis, :]
        
        self.chisq = np.sum(self.chisq_vec) - 2*self.llp
        self.chisq_nd = self.chisq / self.nd
        self.grad = np.sum(grad_vec, axis=1) - 2*self.dllpda
        self.grad_scaled = self.grad * self.param_rescale / self.nd
        self.hessian = np.sum(hess_vec, axis=2) - 2*self.d2llpda2
        self.hessian_scaled = self.hessian * \
            (self.param_rescale * self.param_rescale[:, np.newaxis]) / self.nd

class GHFfitfun:
    ''' Parent class for GHF_xxxx classes.
        '''
    numparams = 0 # class attribute, value varies by child class
    def __init__(self, mask=None, fixedparams = None):
        ''' Inherited by most GHF_xxxx classes.  All GHF_xxxx classes inherit
            the ability to mask some set of parameters, fixing them rather
            than allowing them to float.  By default all parameters float.  If
            mask is supplied, it should be a boolean ndarray with shape
            (numparams, ) -- 'False' elements in the mask will be fixed.
            Fixed elements will be fixed at values specified in fixedparams.
            fixedparams may have shape (numparams, ), or may have a length
            equal to the np.sum(~mask).
            '''
        if mask is None:
            self.parammask  = np.ones((self.numparams), dtype=np.bool)
        else:
            self.parammask = mask.copy()
        self.numparams_eff = np.sum(self.parammask)
        self.params = np.zeros((self.numparams), dtype=np.float64) + np.nan
        if np.any(~self.parammask):
            if len(fixedparams) == self.numparams:
                self.params[~self.parammask] = fixedparams[~self.parammask]
            else:
                self.params[~self.parammask] = fixedparams
        self.binedges = np.zeros((2), dtype=np.float64)
        self.v = np.zeros((self.binedges.size-1), dtype=np.float64)
        self.dvda = np.zeros((self.numparams_eff, self.binedges.size-1),
                              dtype=np.float64)
        self.d2vda2 = np.zeros((self.numparams_eff, self.numparams_eff,
                                self.binedges.size-1),
                               dtype=np.float64)
    
    def update_binedges(self, binedges):
        ''' Updates the binedges and allocates the output arrays accordingly.
            Also runs the calculation if valid values for params exist.
            '''
        self.binedges = binedges.copy()
        self.v = np.zeros((binedges.size-1), dtype=np.float64)
        self.dvda = np.zeros((self.numparams_eff, binedges.size-1),
                              dtype=np.float64)
        self.d2vda2 = np.zeros((self.numparams_eff, self.numparams_eff,
                                binedges.size-1),
                               dtype=np.float64)
        if not np.any(np.isnan(self.params)):
            self.calculate()        

    def update_params(self, newparams):
        ''' Updates the current parameters and runs the calculation --- unless
            the new parameters are identical to the old ones.
            '''
        if not np.all(self.params[self.parammask]==newparams):
            self.params[self.parammask] = newparams
            self.calculate()
    
    def calculate(self):
        ''' This must be filled in for child classes.  This method should
            update v, dvda, and d2vda2 based on values in binedges and params.
            '''
        pass
    
    def cdf(self, x):
        ''' This is a helper function for the default pdf method.
            '''
        oldbinedges = self.binedges.copy()
        self.update_binedges(x)
        result = np.cumsum(self.v.copy())
        self.update_binedges(oldbinedges)
        return result
    
    def pdf(self, x):
        ''' It can be useful to evaluate the 'pdf' (normalized by number of
            entries expected in the histogram) at given points, particularly
            for plotting purposes.  This can be updated in child classes, but
            if not, this will use the cdf to estimate the pdf.
            '''
        x_edges = np.zeros((x.size+1), dtype=np.float64)
        x_edges[1:-1] = 0.5 * np.diff(x) + x[:-1]
        x_edges[0] = 2*x[1]-x[2]
        x_edges[-1] = 2*x[-2] - x[-3]
        result = np.diff(self.cdf(x_edges)) / np.diff(x_edges)
        return result

class GHF_exp(GHFfitfun):
    ''' Exponential distribution.  Defined as:
        a = params[0]
        x0 = params[1]
        P(x) = a * exp(x/x0)
        '''        
    numparams = 2
    def calculate(self):
        a = self.params[0]
        x0 = self.params[1]
        int_y = a * x0 * np.exp(self.binedges/x0)
        self.v = np.diff(int_y)
        self.dvda[0] = self.v / a
        self.dvda[1] = np.diff(int_y * ((1/x0) - (self.binedges/(x0**2))))
        self.d2vda2[1,1] = np.diff(int_y * (self.binedges**2)/(x0**4))
        self.d2vda2[0,1] = self.dvda[1] / a
        self.d2vda2[1,0] = self.d2vda2[0,1]
    def pdf(self, x):
        a = self.params[0]
        x0 = self.params[1]
        result = a * np.exp(x/x0)
        return result

class GHF_flat(GHFfitfun):
    ''' Flat (constant) distribution.  Defined as:
        P(x) = params[0]
        '''
    numparams = 1
    def calculate(self):
        self.dvda[0] = np.diff(self.binedges)
        self.v = self.dvda[0] * self.params[0]
    def pdf(self, x):
        result = self.params[0] + np.zeros(x.shape, dtype=np.float64)
        return result