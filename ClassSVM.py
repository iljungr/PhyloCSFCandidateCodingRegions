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
Support Vector Machine using R package e1071.
Before creating a ClassSVM, the e1071 R package must be installed at least once on the
    machine being used. This can be done by calling install_e1071 from an interactive
    session.
"""
from __future__ import division
from __future__ import print_function
import sys

# Don't import rpy2 here, so scripts importing this file needn't import rpy2 unless they
# actually invoke these utilities. Instead, import rpy2 in the actual methods.

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

class ClassSVM(object) :
    """
    Support Vector Machine.
    To use, create an SVM class instance using training points and then get SVM scores
        by calling the SVM class instance with a list of test points.
    trainingPoints:  iterable of one or more k-dimensional points (all with the same k)
    trainingClasses: the corresponding classes (integers starting at 0)
    testPoints:      iterable of one or more k-dimensional points (same k as training)
    kArgs are additional keyword arguments to be passed to R's svm function. See
        https://cran.r-project.org/web/packages/e1071/e1071.pdf for arguments and values.
        Note that default radial basis kernel is 'radial', not 'radial basis'.
        Example: ClassSVM(trP, trC, kernel = 'radial', cost = 2, probability = False)
    Note that even when probability is set to True (the default), the "probabilities" it
        returns are sometimes a bit less than 0 or greater than 1.
    """
    def __init__(self, trainingPoints, trainingClasses, **kArgs) :
        """
        Train an SVM using the training points with default parameters unless overridden
            by kArgs.
        """
        import rpy2.robjects as robjects
        import rpy2.robjects.packages as rpackages
        rpackages.importr('e1071')
        self.svmFunc = robjects.r['svm'](r_matrix_from_points(trainingPoints),
                                         robjects.IntVector(trainingClasses),
                                         **kArgs)
    def __call__(self, testPoints) :
        """
        Return the list of class "probabilities" of the test points.
        """
        import rpy2.robjects as robjects
        return list(robjects.r['predict'](self.svmFunc, r_matrix_from_points(testPoints)))

def r_matrix_from_points(points) :
    """
        Given a list (or other iterable) of one or more k-dimensional points (all the same k)
        return a representation of it as an rpy2 robject of type matrix
        """
    import rpy2.robjects as robjects
    combinedPoints = [x for point in points for x in point]
    return robjects.r.matrix(robjects.FloatVector(combinedPoints), nrow = len(points),
                             byrow = True)

def install_e1071() :
    # This only needs to be called once for each package on each machine.
    # Do it from an interactive session, because it might ask questions.
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects import StrVector
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    utils.install_packages(StrVector(['e1071']))
