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

#include <vector>

#include "ess.hh"
#include "quality_box.hh"

void BoxQualityFunction::create_integral_matrices(const std::vector<double> &raw_matrix) {
    pos_matrix.clear();
    pos_matrix.resize( raw_matrix.size() );
    neg_matrix.clear();
    neg_matrix.resize( raw_matrix.size() );

    // split the weight matrix into positive and negative entries
    for (unsigned int i=0; i < raw_matrix.size(); i++) {
        double val = raw_matrix[i];
        if (val > 0.)
            pos_matrix[i] = val;
        else
            neg_matrix[i] = val;
    }

    // calculate integral image verically
    for (int j=1; j < height; j++) {
        for (int i=1; i < width; i++) {
            pos_matrix[off(i,j)] += pos_matrix[off(i,j-1)];
            neg_matrix[off(i,j)] += neg_matrix[off(i,j-1)];
        }
    }
    // calculate integral image horizontally
    for (int j=1; j<height; j++) {
        for (int i=1; i<width; i++) {
            pos_matrix[off(i,j)] += pos_matrix[off(i-1,j)];
            neg_matrix[off(i,j)] += neg_matrix[off(i-1,j)];
        }
    }
    return;
}

void BoxQualityFunction::setup(int argnumpoints, int argwidth, int argheight, 
                               double* argxpos, double* argypos, double* argclst, 
                               void* argdata) {
    width = argwidth;
    height = argheight;
    
    // transform (x,y,c),weight into integral image representation
    std::vector<double> raw_matrix(argwidth*argheight);

    // for sum-of-scores, the data is a vector of cluster weights
    const double* argweight = reinterpret_cast<double*>(argdata);

    // we pad +1 so we can avoid boundary checks later
    for (int k=0; k<argnumpoints; k++) {
        const int x = static_cast<int>(argxpos[k])+1;
        const int y = static_cast<int>(argypos[k])+1;
        const int c = static_cast<int>(argclst[k]);
        raw_matrix[off(x,y)] += argweight[c];
    }
    create_integral_matrices(raw_matrix);
    return;
}

void BoxQualityFunction::cleanup() {
    return;
}

double BoxQualityFunction::upper_bound(const sstate* state) const {
    return quality_upper_single(state);
}

