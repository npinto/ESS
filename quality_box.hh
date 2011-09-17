/********************************************************
 *                                                      *
 *  Efficient Subwindow Search (ESS) implemented in C++ *
 *  bound the sum-of-entries in a box (e.g. linear SVM) *
 *                                                      *
 *   Copyright 2006-2008 Christoph Lampert              *
 *   Contact: <mail@christoph-lampert.org>              *
 *                                                      *
 *  Licensed under the Apache License, Version 2.0 (the *
 *  "License"); you may not use this file except in     *
 *  compliance with the License. You may obtain a copy  *
 *  of the License at                                   *
 *                                                      *
 *     http://www.apache.org/licenses/LICENSE-2.0       *
 *                                                      *
 *  Unless required by applicable law or agreed to in   *
 *  writing, software distributed under the License is  * 
 *  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES *
 *  OR CONDITIONS OF ANY KIND, either express or        *
 *  implied. See the License for the specific language  *
 *  governing permissions and limitations under the     *
 *  License.                                            *
 *                                                      *
 ********************************************************/

#ifndef _QUALITY_BOX_H
#define _QUALITY_BOX_H

#include <vector>

#include "ess.hh"
#include "quality_function.hh"

class BoxQualityFunction : public QualityFunction {

    private:
        int width,height;
        std::vector<double> pos_matrix;
        std::vector<double> neg_matrix;

        // convert (x,y) into 1d index
        inline unsigned int off(unsigned int x, unsigned int y) const {
            return y*width+x;
        }
        
        // calculate score of a box from integral image
        double rect_val(unsigned int xl, unsigned int yl, 
                        unsigned int xh, unsigned int yh,
                        const std::vector<double> &matrix) const {
            if ( (xl > xh) || (yl > yh)) return 0.;

            const double val = matrix[off(xh,yh)] - matrix[off(xh,yl-1)]
                               - matrix[off(xl-1,yh)] + matrix[off(xl-1,yl-1)];
            return val;
        }

        // calculate upper bound for one set of rectangles
        double quality_upper_single(const sstate* s) const {
            const double fplus = rect_val(s->low[0], s->low[1], s->high[2], s->high[3], pos_matrix);
            const double fminus = rect_val(s->high[0], s->high[1], s->low[2], s->low[3], neg_matrix);
            return fplus+fminus;
        }


        // create separate integral images for positive and negative part
        // of the original weight matrix
        void create_integral_matrices(const std::vector<double> &raw_matrix);

    public:
        void setup(int argnumpoints, int argwidth, int argheight, 
                   double* argxpos, double* argypos, double* argclst, 
                   void* argdata);

        void cleanup();

        double upper_bound(const sstate* state) const;
};
#endif
