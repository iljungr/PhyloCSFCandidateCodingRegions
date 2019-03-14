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
MinimizeUtil.py
Utility for minimizing a function from R^n to R.
"""
from __future__ import division, print_function
import sys, itertools, numpy
from math import sqrt
from pprint import pprint

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

def minimize(f, guess, xscales, relxtol, maxNumSteps = None, dump = False) :
    """
    Find the point that minimizes f:R^n->R given a starting guess.
    Return x, f(x)
    xscales is a vector specifying the rough "scale" of each coordinate. The initial
      simplex is constructed using these scales as the sizes of the legs. It is considered
      converged roughly when each coordinate of the result is off by no more than
      relxtol times the corresponding coordinate of xscale.
    """
    n = len(guess)
    if maxNumSteps == None :
        maxNumSteps = 30 * n * n
    initialSimplex = [list(guess)]
    for ii in range(n) :
        vertex = list(guess) # Copy, and make sure it is mutable.
        vertex[ii] += xscales[ii]
        initialSimplex.append(vertex)

    def f_exception(x) :
        try :
            return f(x)
        except :
             print('Failed to evaluate', x, file = sys.stderr)
             raise

    return Nelder_Mead(f_exception, initialSimplex,
                       [relxtol * xscale for xscale in xscales],
                       maxNumSteps, dump)

class NoConvergenceError(Exception) : pass

def Nelder_Mead(f, initialSimplex, xtolVec = None, maxNumSteps = None, 
                dump = False,  # Debugging option
                ) :
    """
    Perform the Nelder-Mead minimization method to find a minimum of f:R^n->R.
    initialSimplex is a sequence of n + 1 n-dimensional points.
    xtolVec is a sequence of tolerances, one for each dimension. Stop when the diameter
        of the simplex in every dimension is less than the corresponding tolerance,
        provided the last step was not an expansion (o.w. it could converge on first step
        if input is very small) or a reduction (because almost all new).
    If number of step exceeds maxNumSteps, raise NoConvergenceError 
    Return x, f(x) at the minimum.
    """
    n = len(initialSimplex) - 1
    
    if xtolVec == None :
        xtolVec = [sqrt(1e-9)] * n
    if maxNumSteps == None :
        maxNumSteps = 30 * n * n
    
    assert all(len(x) == n for x in initialSimplex)
    simplex = [(numpy.array(map(float, x)), f(x)) for x in initialSimplex]
    simplex.sort(key = lambda pair : pair[1]) 
    
    prevStepExpansionOrReduction = True # Disallow convergence on first step
    
    for count in itertools.count() :

        if dump :
            print('Simplex:', file = sys.stderr)
            pprint([(map(lambda t : '%.5g' % t, list(x)), '%.5g' % fx)
                    for x, fx in simplex],
                   sys.stderr)
        
        # Test stopping condition
        if not prevStepExpansionOrReduction :
            for ii in range(n) :
                minX = min(x[ii] for x, fx in simplex)
                maxX = max(x[ii] for x, fx in simplex)
                if maxX - minX >= xtolVec[ii] :
                    break
            else :
                if dump :
                    print('Converged in %d steps.' % count, file = sys.stderr)
            
                return simplex[0] # x, f(x)
        if count >= maxNumSteps :
            raise NoConvergenceError
            
        stepType = simplex_one_step(f, simplex)
        
        prevStepExpansionOrReduction = (stepType in ['Expansion', 'Reduction'])
        
        if dump :
            print('Choosing', stepType, file = sys.stderr)
         
def simplex_one_step(f, simplex) :
    # Perform one step of Nelder-Mead minimization on the simplex [(x, f(x)), ...]
    #    (changing it in place).
    # Return string stating which kind of step was taken.
    # simplex is assumed to be sorted by increasing f(x) on input, and will be on return.
    n = len(simplex) - 1
    centroid = sum([x for x, fx in simplex[:-1]], numpy.zeros(n)) / n
    reflection = centroid + (centroid - simplex[n][0])
    freflection = f(reflection)
    if simplex[0][1] <= freflection < simplex[n - 1][1] :
        stepType = 'Reflection'
        simplex[n] = (reflection, freflection)
    elif freflection < simplex[0][1] :
        expansion = centroid + 2 * (centroid - simplex[n][0])
        fexpansion = f(expansion)
        if fexpansion < freflection :
            stepType = 'Expansion'            
            simplex[n] = (expansion, fexpansion)
        else :
            stepType = 'Reflection'            
            simplex[n] = (reflection, freflection)
    else :
        contraction = centroid - .5 * (centroid - simplex[n][0])
        fcontraction = f(contraction)
        if fcontraction < simplex[n][1] :
            stepType = 'Contraction'            
            simplex[n] = (contraction, fcontraction)
        else :
            stepType = 'Reduction'            
            for ii in range(1, n + 1) :
                newx = simplex[0][0] + 0.5 * (simplex[ii][0] - simplex[0][0])
                simplex[ii] = (newx, f(newx))
    simplex.sort(key = lambda pair : pair[1]) # Could improve by inserting at right place...
    return stepType
