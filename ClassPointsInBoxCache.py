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
from __future__ import division
import sys, itertools

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

## Interface:
class PointsInBoxCache :
    # This class is a cache allowing fast determination of which of a set of 2D points lie in a box.
    # Only member variable: quadtree (for internal use)
    def __init__(self, points = []) :
        # points is an iterator of objects to which [0] and [1] can be applied to get x and y.
        # For example, they could be 2-tuples (x, y); but they could also be more complex, 
        #     such as lists [x, y, data].
        # The points will be referenced but not copied,
        # Duplicate points will NOT be ignored.
        init_cache(self, points)
        
    def add(self, point) :
        # Note that the cache will work more efficiently if all points are added in __init__.
        pibc_add(self, point)
        
    def remove(self, point) :
        # Remove the first point found that matches point (using "==") or raise ValueError.
        pibc_remove(self, point)
        
    def points_in_box(self, xmin, xmax, ymin, ymax) :
        # Return an iterator through all of the initial points
        #   that satisfy xmin <= point[0] <= xmax and ymin <= point[1] <= ymax.
        # Note that if a pair appeared multiple times in the original list, it will
        #   be included multiple times in the result.
        return points_in_box(self, xmin, xmax, ymin, ymax)


## Implementation:

QUAD_TREE_SPLIT_SIZE = 50 # Determined empirically (though you'd want it bigger if most 
                          #   points_in_box queries find huge number of points)
class QuadTree :
    """
    points: list of <= QUAD_TREE_SPLIT_SIZE points 
    mid: [midx, midy]
    children: [[LeftBottomQuadTree, LeftTopQuadTree], [RightBottomQuadTree, RightTopQuadTree]]
    If self contains no points, everything is None.
    If self contains 1 to QUAD_TREE_SPLIT_SIZE points, then children and mid are None.
        (Actually, node with 1 to QUAD_TREE_SPLIT_SIZE can have children after remove is called.)
    If self contains more than QUAD_TREE_SPLIT_SIZE points, then points is None.
    All points in children[0][0] have point[0] <= midx and point[1] <  midy.
    All points in children[0][1] have point[0] <= midx and point[1] >= midy.
    All points in children[1][0] have point[0] >= midx and point[1] <  midy.
    All points in children[1][1] have point[0] >= midx and point[1] >= midy.
    (Points with point[0] == midx can be either left or right, but point[1] is different)
    ...
    """
    def __init__(self, points) :
        # points are any objects for which [0] and [1] yield numbers. Needn't be unique.
        self.points = self.mid = self.children = None
        if len(points) == 0 :
            return
        self.points = points # Reuses list! Should be OK in current implementation.
        self.split_if_too_big()
            
    def split_if_too_big(self) :
        points = self.points
        if len(points) <= QUAD_TREE_SPLIT_SIZE :
            return
        centerInd = len(points) // 2
        points.sort(key = lambda pt : (pt[0], pt[1])) # Don't peek at anything but [0] and [1]
        midx = (points[centerInd - 1][0] + points[centerInd][0])  / 2
        # Instead of just dividing the x points based on the mid value, 
        #   take half of the sorted list. This  guarantees that number of points at each level
        #   decreases exponentially even if there are lots of equalities.
        leftPoints = points[:centerInd]
        rightPoints = points[centerInd:]
        # Can't do same for y, but don't need to.
        points.sort(key = lambda pt : (pt[1], pt[0])) # In order to find median
        midy = (points[centerInd - 1][1] + points[centerInd][1])  / 2
        self.children = [
                [  QuadTree([pt for pt in leftPoints  if pt[1] <  midy]),
                   QuadTree([pt for pt in leftPoints  if pt[1] >= midy])  ],
                [  QuadTree([pt for pt in rightPoints if pt[1] <  midy]),
                   QuadTree([pt for pt in rightPoints if pt[1] >= midy])  ],
            ]
        self.mid = [midx, midy]
        self.points = None
        
    def add(self, point) :
        if self.mid != None :
            # If point[0] == midx, it could go in 2 places; ideally, I'd keep counts
            # for each child and put it in the smaller one...
            self.children[point[0] > self.mid[0]][point[1] >= self.mid[1]].add(point)
        else :
            if self.points == None :
                self.points = []
            self.points.append(point)
            self.split_if_too_big()
    
    def remove(self, point) :
        # Remove the first point found that matches point (using "==") or raise ValueError.
        # Clean up nodes with all empty children, but don't rebalance tree, and don't
        #     move  points from children back to self when count drops below QUAD_TREE_SPLIT_SIZE.
        if self.points != None :
            del self.points[self.points.index(point)] # Rasies ValueError if not found.
            if len(self.points) == 0 :
                self.points = None
            return
        if self.children == None :
            raise ValueError
        xinds = [0] if point[0] < self.mid[0] else [1] if point[0] > self.mid[0] else [0, 1]
        found = False
        for xind in xinds :
            child = self.children[xind][point[1] >= self.mid[1]]
            try :
                child.remove(point)
                found = True
            except ValueError :
                pass
        if not found :
            raise ValueError
        if all(child.points == None and child.children == None 
               for child in itertools.chain.from_iterable(self.children)) :
            self.children = self.mid = None
        
    def iterPointsInBox(self, xmin, xmax, ymin, ymax) :
        # Iterate through all points satisfying xmin <= x <= xmax and ymin <= y <= ymax
        if self.points != None :
            for point in self.points :
                if xmin <= point[0] <= xmax and ymin <= point[1] <= ymax :
                    yield point
        if self.children == None :
            return
        goodxinds, goodyinds = [], []
        if self.mid[0] >= xmin :
            goodxinds.append(0)
        if self.mid[0] <= xmax :
            goodxinds.append(1)
        if self.mid[1] >= ymin :
            goodyinds.append(0)
        if self.mid[1] <= ymax :
            goodyinds.append(1)
        for xind in goodxinds :
            for yind in goodyinds :
                for point in self.children[xind][yind].iterPointsInBox(xmin, xmax, ymin, ymax) :
                    yield point

def init_cache(cache, points) :
    cache.quadtree = QuadTree(list(points))

def pibc_add(cache, point) :
    cache.quadtree.add(point)

def pibc_remove(cache, point) :
    cache.quadtree.remove(point)
    
def points_in_box(cache, xmin, xmax, ymin, ymax) :
    return cache.quadtree.iterPointsInBox(xmin, xmax,  ymin, ymax)
    
# For debugging:
def print_qtree(qtree, indent = 0) :
    tabs = '\t' * indent
    print tabs + 'points =', qtree.points
    print tabs + 'mid = ', qtree.mid
    if qtree.children != None :
        print tabs + 'LeftBottom:'
        print_qtree(qtree.children[0][0], indent + 1)
        print tabs + 'LeftTop:'
        print_qtree(qtree.children[0][1], indent + 1)
        print tabs + 'RightBottom:'
        print_qtree(qtree.children[1][0], indent + 1)
        print tabs + 'RightTop:'
        print_qtree(qtree.children[1][1], indent + 1)
        
def qtree_stats(qtree) :
    # Return number of nodes, number of points
    if qtree == None :
        return 0, 0
    if qtree.points != None :
        assert(qtree.children == None and qtree.mid == None)
        return 1, len(qtree.points)
    if qtree.children != None :
        assert(qtree.mid != None)
        numNodes = 1
        numPoints = 0
        for child in itertools.chain.from_iterable(qtree.children) :
            addNodes, addPoints = qtree_stats(child)
            numNodes += addNodes
            numPoints += addPoints
        return numNodes, numPoints
    return 1, 0
