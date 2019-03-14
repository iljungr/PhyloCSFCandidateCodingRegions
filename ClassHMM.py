#!/usr/bin/env python
# Copyright 2015 by Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Utilities for computing best path and state probabilities for an HMM (Hidden Markov Model)
    with known initial, transition, and emission probabilities.
HMM can have continuous emissions (class HMM) or discrete emissions (class FiniteEmissionHMM).
"""
from __future__ import division
import numpy, array, sys
from math import log

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

def almost1(x) :
    return 1 - 1e-9 < x < 1 + 1e-9

class HMM(object) :
    def __init__(self, numStates, initProbs, transProbs, emitProb) :
        """
        numStates = Number of hidden states. Hidden states are named 0, 1, ... numStates-1.
        initProbs = Prob(starting in state i)
        transProbs[i][j] = Prob(state i -> state j)
        emitProb(s, e) = Prob(state s emits e)
        Note that for most purposes it is OK if for a given emission emitProbs is multiplied
            by a (positive) factor that is constant across states. This will not affect best_path,
            or state_probabilities, and it will change log_prob_obs by adding logs of these factors 
            for all positions. Similarly, emitProbs can be a probability density rather than a 
            probability.
        """
        self.numStates = numStates
        self.initProbs = initProbs
        self.transProbs = transProbs
        self.emitProb = emitProb
        assert(len(initProbs) == numStates)
        assert(len(transProbs) == numStates)
        assert(all(len(row) == numStates for row in transProbs))
        assert(almost1(sum(initProbs)))
        assert(all(almost1(sum(row)) for row in transProbs))
        assert(all(p >= 0 for p in initProbs))
        assert(all(p >= 0 for row in transProbs for p in row))

    def log_prob_obs(self, observation) :
        "observation is a sequence of emissions. Return log prob of seeing it under this hmm."
        observation = iter(observation) # Handle list, iterator, etc.
        initObs = observation.next()
        logFactor = 0
        prevProbs = [initProb * self.emitProb(state, initObs) 
                     for state, initProb in enumerate(self.initProbs)]
        for curObs in observation :
            # exp(logFactor) * curProbs are the probabilities of being in a particular state
            # at this point and having seen all observations up to and including this point
            curProbs = [sum(prevProbs[prevState] * self.transProbs[prevState][curState] 
                            for prevState in range(self.numStates)
                           ) * self.emitProb(curState, curObs) 
                        for curState in range(self.numStates)]
            maxProb = max(curProbs)
            logFactor += log(maxProb)
            curProbs = [p / maxProb for p in curProbs]
            prevProbs = curProbs
        return logFactor + log(sum(prevProbs))

    def best_path(self, observation) :
        """Return the most likely sequence of hidden states given the observation,
           using the Viterbi algorithm."""
        hiddenStates = range(self.numStates)
        observation = iter(observation) # Handle list, iterator, etc.
        initObs = observation.next()
        
        # bestStates[j][k] is the hidden state at position j that maximized the
        #    probability of the observation, given hidden state k at position j + 1
        bestStates = []
        
        # prevBestProbs[k] is the probability of the best path from the start to the
        #    previous position with state k at the previous position, 
        #    times a constant factor independent of k
        prevBestProbs = [initProb * self.emitProb(state, initObs) 
                         for state, initProb in enumerate(self.initProbs)]
                         
        for curObs in observation :
            bestStates.append(array.array('i'))
            curBestProbs = []
            for curState in hiddenStates :
                newProb, bestPrevState = max((prevBestProbs[prevState] * 
                                              self.transProbs[prevState][curState],
                                              prevState)
                                             for prevState in hiddenStates)
                newProb *= self.emitProb(curState, curObs)
                bestStates[-1].append(bestPrevState)
                curBestProbs.append(newProb)
            # Scale to avoid floating underflow
            maxProb = max(curBestProbs)
            curBestProbs = [p / maxProb for p in curBestProbs]
            prevBestProbs = curBestProbs
        result = []
        result.append(indexmax(prevBestProbs))
        for bestStatesForPos in bestStates[::-1] :
            result.append(bestStatesForPos[result[-1]])
        return result[::-1]
        
    def state_probabilities(self, observationSeq) :
        """ Return result, where result[pos][s] = probability position pos is in state s.
            observationSeq must be a sequence, not an iterator; we will traverse both
            directions."""
        numObs = len(observationSeq)
        numStates = self.numStates
        transProbs = self.transProbs
        emitProb = self.emitProb
        states = range(numStates)
        
        # forward[pos][s] = (prob of being in state s at position pos and emitting
        #     observationSeq[0 ... pos]) * ArbitraryScaleFactor(pos)
        # ArbitraryScaleFactor is a scale factor to avoid floating overflow. As long as we
        #     apply the same scale factor to forward[pos][s] for all s, it doesn't change
        #     the outcome.
        forward = numpy.zeros([numObs, numStates])
        for state in states :
            forward[0][state] = self.initProbs[state] * emitProb(state, observationSeq[0])
        for pos in xrange(1, len(observationSeq)) :
            obs = observationSeq[pos]
            for state in states :
                forward[pos][state] = \
                    sum(forward[pos - 1][prevState] * transProbs[prevState][state]
                        for prevState in states) * emitProb(state, obs)
            maxf = max(forward[pos])
            for state in states :
                forward[pos][state] /= maxf
        
        # backward[pos][s] = (prob of emitting observationSeq[pos + 1 ... -1] given
        #     being in state s at position pos) * ArbitraryScaleFactor(pos)
        # ArbitraryScaleFactor is a scale factor to avoid floating overflow. As long as we
        #     apply the same scale factor to backward[pos][s] for all s, it doesn't change
        #     the outcome.
        backward = numpy.zeros([numObs, numStates])
        for state in states :
            backward[-1][state] = 1
        for pos in xrange(len(observationSeq) - 2, -1, -1) :
            obs = observationSeq[pos + 1]
            for state in states :
                backward[pos][state] = \
                    sum(transProbs[state][nextState] * emitProb(nextState, obs) *
                        backward[pos + 1][nextState]
                        for nextState in states)
            maxb = max(backward[pos])
            for state in states :
                backward[pos][state] /= maxb
                
        result = numpy.zeros([numObs, numStates])
        prods = numpy.zeros(numStates) # Reuse to avoid memory allocation
        for pos in xrange(len(observationSeq)) :
            for state in states :
                prods[state] = forward[pos][state] * backward[pos][state]
            total = sum(prods)
            for state in states :
                result[pos][state] = prods[state] / total

        return result
                        
class FiniteEmissionHMM(HMM) :
    def __init__(self, numStates, initProbs, transProbs, emitProbs) :
        """ Same as HMM, except that the emission probabilities are specified by a matrix
            emitProbs[i][j] = P(state i emits j)"""
        assert(len(set(len(row) for row in emitProbs)) == 1) # All rows have same length
        # assert(all(almost1(sum(row)) for row in emitProbs)) # Commented out to allow scaling
        assert(all(p >= 0 for row in emitProbs for p in row))
        assert(len(emitProbs) == numStates)
        self.emitProbs = emitProbs # May be useful for debugging
        super(FiniteEmissionHMM, self).__init__(numStates, initProbs, transProbs,
            (lambda state, emits : self.emitProbs[state][emits]))

def indexmax(iter) :
    # Return the index of the maximum element in seq
    bestInd = None
    for ind, elt in enumerate(iter) :
        if bestInd == None or elt > maxElt :
            maxElt = elt
            bestInd = ind
    return bestInd     
