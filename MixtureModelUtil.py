#!/usr/bin/env python
# Copyright 2019 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
MixtureModelUtil.py.
Utility for optimizing the parameters of a mixture model.
"""
from __future__ import division, print_function
from math import log, exp, sqrt, pi
import sys
from MinimizeUtil import minimize

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

class ParametrizedLogDistribution(object) :
    # Trivial base class defining a family of distributions controlled by a sequence of parameters.
    # Real class must define its own function and gradient.
    def log_density_and_grad(self, point, params) :
        # Return the log_density of the distribution at point for the given parameters.
        # Point could be a number, a sequence, or whatever.
        return 0.0

def infer_mixture(points, paramDists, paramGuesses, classPriorGuesses,
                  numSteps = 10, relxtol = 0.001, dump = False) :
    """
    Use Expectation Maximization to find the parameter values and priors for a set of 
        distributions to define a mixture model that maximizes the expectation of sampling the 
        points.
    Inputs:
      points: A sequence of points. Each can be represnted as a list, a number, or anything else
        that can be passed as the first argument to the ParametrizedDistribution methods.
      paramDists: a sequence of ParametrizedLogDistribution classes defining the families of distributions.
        The number of parameters for each one is defined by the length of the corresponding
        sequence in parameterGuesses. 
      paramGuesses: a set of sequences, each of which is the initial guess for the
        parameters of the corresponding function.
      classPriorGuesses: initial guesses for the class priors. Needn't sum to 1 -- they'll be
        normalized.
      numSteps: number of EM steps
      relxtol: relative tolerance of the parameter values, used for termination criterion
        of minimize
    
    Return:
      parameters: sequence of sequences, same size as parameter guesses
      classPriors: sequence of numbers, one for each distribution
      classProbs: [[prob point0 in class0, prob point0 in class1, ...], ...]
    """
    numPoints = len(points)
    numDist = len(paramDists)
    assert(numDist == len(paramGuesses) == len(classPriorGuesses))
    paramSeqs = [list(guesses) for guesses in paramGuesses]
    classPriors = list(classPriorGuesses)
    classProbs = [[0.0] * numDist for ii in range(numPoints)] # To avoid repeated reallocation
    
    for iteration in range(numSteps) : # Until a better stopping condition is found...
        if not dump :
            print('EM iteration', iteration, numSteps, file = sys.stderr)
        # E step: find class probabilities
        for ii in range(numPoints) :
            likelihoods = [classPriors[jj] * 
                           exp(paramDists[jj].log_density_and_grad(points[ii], paramSeqs[jj]))
                           for jj in range(numDist)]
            total = sum(likelihoods)
            for jj in range(numDist) :
                classProbs[ii][jj] = likelihoods[jj] / total if total != 0 else 1 / numDist
                
        # M step: find best parameter values and class priors given class probabilities
          
        """
        The log likelihood is a sum in which the terms depending on classPriors separate out and
           can be minimized analytically. The result is what you'd expect: the average classProbs.
        Note that the latent variable is the discrete class _membership_, not the class 
            probabilities. Thus, we aren't figuring out the likelihood of having particular class 
            probabilities, we are using the class probabilities as weights on the likelihood of 
            particular class membership. Think of the class probabilities as defining
            fractional points, with each point in one discrete class or another.
        """
        classPriors[:] = [sum(classProbs[ii][jj] for ii in range(numPoints)) / numPoints
                          for jj in range(numDist)]

        # Find best params for each distribution (they can be done separately because
        # log likelihood is a sum of terms each of which depends on only one of them).
        for jj in range(numDist) :
            if len(paramSeqs[jj]) == 0 :
                continue # This distribution is not parametrized. Nothing to minimize.
            def fjj(paramsJJ) : # - log likelihood
                fOfParams = 0
                for ii in range(numPoints) :
                    logLikelihood = paramDists[jj].log_density_and_grad(points[ii], paramsJJ)
                    fOfParams -= classProbs[ii][jj] * logLikelihood
                return fOfParams
        
            paramSeqs[jj] = list(minimize(fjj, paramSeqs[jj],
                                          xscales = [0.1] * len(paramSeqs[jj]),
                                          relxtol = relxtol, dump=False)[0])
        if dump :
            print(iteration, paramSeqs, classPriors, file = sys.stderr)
         
    return paramSeqs, classPriors, classProbs
