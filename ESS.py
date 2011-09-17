#!/usr/bin/python

# ********************************************************
# *                                                      *
# * Efficient Subwindow Search (ESS) from Python         *
# *                                                      *
# * Wrapper to call ESS from Python (using ctypeslib)    *
# *                                                      *
# *   Copyright 2006-2008 Christoph Lampert              *
# *   Contact: <mail@christoph-lampert.org>              *
# ********************************************************

from ctypes import Structure,c_int,c_double
from numpy.ctypeslib import load_library,ndpointer

class Box_struct(Structure):
        """Structure to hold left,top,right,bottom and score of a box instance.
           The fields have to coincide with the C-version in pyramid_search.h
        """
        _fields_ = [("left", c_int), 
                    ("top", c_int), 
                    ("right", c_int), 
                    ("bottom", c_int), 
                    ("score", c_double) ]


def subwindow_search_pyramid(numpoints, width, height, xpos, ypos, clstid, numbins, numlevels, weights):
    """Subwindow search for best box with bag-of-words histogram with a spatial 
       pyramid kernel."""
    pyramidlib = load_library("libess.so",".")

    pyramidlib.pyramid_search.restype = Box_struct
    pyramidlib.pyramid_search.argtypes = [c_int,c_int,c_int,        
            ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS'),
            ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS'),
            ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS'),
            c_int, c_int,                                           
            ndpointer(dtype=c_double, ndim=1, flags='C_CONTIGUOUS')]

    box = pyramidlib.pyramid_search(numpoints, width, height, 
                      xpos, ypos, clstid, numbins, numlevels, weights)
    return box


def subwindow_search(numpoints, width, height, xpos, ypos, clstid, weights):
    """Subwindow search for best box with linear bag-of-words kernel."""
    return subwindow_search_pyramid(numpoints, width, height, xpos, ypos, \
                                    clstid, max(clstid)+1, 1, weights)


# Example of usage: load x,y,clst and weight files and search for best box.
if __name__ == "__main__":
    import sys
    import pylab
    from numpy import array
    
    try:
        xyc = pylab.load(sys.argv[1])
        x=array(xyc[:,0], c_double).copy() # force to be contiguous
        y=array(xyc[:,1], c_double).copy()
        c=array(xyc[:,2], c_double).copy()
        w = pylab.load(sys.argv[2])
        numpoints = len(x)
        width = max(x)+10
        height = max(y)+10
    except IndexError:
        print "Usage: %s featurefile weightfile [number-of-pyramid-levels]\n"
        raise SystemExit
    except IOError:
        print "Can't open input files.\n"
        raise SystemExit

    try:
        numlevels = int(sys.argv[3])
    except IndexError: 
        numlevels = 1
    numbins = w.shape[0]/(numlevels*(numlevels+1)*(2*numlevels+1)/6)

    bestBox = subwindow_search_pyramid(numpoints, int(width), int(height), \
                                       x, y, c, numbins, numlevels, w)

    print "box found: [left: %d, top: %d, right: %d, bottom: %d, score: %f]" \
             % (bestBox.left,bestBox.top,bestBox.right,bestBox.bottom,bestBox.score)

