 ********************************************************
 *                                                      *
 * Efficient Subwindow Search (ESS) implemention in C++ *
 * entry point and commmand line interface              *
 *                                                      *
 * version 1.0, March 27th 2008                         *
 *   first release, except bugfixes soon                *
 *                                                      *
 *   Copyright 2006-2008 Christoph Lampert              *
 *   Contact: <mail@christoph-lampert.org>              *
 *                                                      *
 *  Performs a branch-and-bound search for an arbitrary *
 *  quality function, provided a routine to bound it:   *
 *  see  C. H. Lampert, M. B. Blaschko and T. Hofmann:  *
 *       Beyond Sliding Windows: Object Localization    *
 *       by Efficient Subwindow Search, CVPR (2008)     *
 *                                                      *
 ********************************************************

This is sample implementation of Efficient Subwindow Search,
a method to replace sliding windows approaches to object
localization by a branch-and-bound search. 

The default quality function is to use a spatial pyramid kernel 
with levels of size 1x1, 2x2 ... NxN where N can be chosen at 
runtime. Default is N=1, i.e. bag-of-visual-words model. 



USAGE: 

./ess imagewidth imageheight weight-file data-file

The weight-file is an ASCII file of 1 entry per line. 
It contains the concatenation of the per-cluster weights for 
each pyramid level. The levels are ordered by increasing index.
Within each level, weights are aranged in reading order. All 
weights must be present, as the number of cluster centers 
is calculated by dividing the total number of weights by the
total number of pyramid cells, which is N*(N+1)*(2N+1)/6. 

The data-file is an ASCII File of 3 columns. Each line has the format

X Y CLST

where (X,Y) are the coordinates of a feature in the image and CLST is
an ID attached to it, e.g. it's ID in a codebook representation.
The ID and X,Y are used to look up the correct weight for this feature 
point in a candidate box.

Output is a single line with all detected box in the format

SCORE LEFT TOP RIGHT BOTTOM 

repeating as often as there are boxes.


The main routine 

Box pyramid_search(int argnumpoints, int argwidth, int argheight, 
                   double* argxpos, double* argypos, double* argclst,
                   int argnumclusters, int argnumlevels, double* argweight)

can also be called from an external application. Extracted multiple boxes 
requires a loop, see the main() routine.



EXAMPLE FILES: 

examples/cow.weight is a weightfile with 3000 clusters in 1 level.
examples/car-l1.weight is a weightfile with 1000 clusters in 1 level.
examples/car-l2.weight is a weightfile with 1000 clusters in 2 levels, 
making it 5000 weights in total.

examples/cow.clst is a data-file of 13498 feature points.
examples/cow.clst is a data-file of 13498 feature points.

Arguments can be passed to the executable by environment variables:

numlevels=2 ./ess 151 101 examples/car-l2.weight examples/car.clst

returns the best box for a 2-level pyramid. 

maxresults=3 ./ess 151 101 examples/car-l1.weight examples/car.clst

returns the 3 best car detections instead of just 1. Both arguments
can be combined:

maxresults=2 numlevels=2 ./ess 151 101 examples/car-l2.weight examples/car.clst


Outputs for the examples are:

./ess 368 272 examples/cow.weight examples/cow.clst
# 10978.8779297 89 117 258 209 
	
./ess 151 101 examples/car-l1.weight examples/car.clst
# 1.51183950901 39 56 118 77 0.286267101765 4 61 144 65 0.138369500637 53 88 147 90 

numlevels=2 ./ess 151 101 examples/car-l2.weight examples/car.clst
# 1.69141745567 33 55 115 77 

maxresults=2 numlevels=2 ./ess 151 101 examples/car-l2.weight examples/car.clst
# 1.69141745567 33 55 115 77 0.15714520216 4 60 136 65 



TEST FILES: 

Apart from the examples, there is one synthetic test case so far: 

maxresults=4 ./ess 5 5 examples/test_corners.weight examples/test_corners.clst

should return 

# 1.39999997616 4 4 4 4 1.29999995232 4 0 4 1 1.20000004768 1 4 1 4 1.10000002384 1 1 1 1 

