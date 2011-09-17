/********************************************************
 *                                                      *
 *  Efficient Subwindow Search (ESS) implemented in C++ *
 *  bound for sum of grid cells, e.g. spatial pyramid   *
 *  a grid is really just a collection of boxes with    *
 *  right way to access them and add up their scores    *
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
#include "quality_pyramid.hh"

sstate PyramidQualityFunction::rel_to_abs_coordinate(const Cell &subcoordinate, 
                                                     const sstate* fullstate) const {
    sstate substate;
    substate.upper = fullstate->upper;
    substate.low[0] =  static_cast<short>((1-subcoordinate.left)  *fullstate->low[0]  + subcoordinate.left  *fullstate->low[2] );
    substate.high[0] = static_cast<short>((1-subcoordinate.left)  *fullstate->high[0] + subcoordinate.left  *fullstate->high[2]);
    substate.low[2] =  static_cast<short>((1-subcoordinate.right) *fullstate->low[0]  + subcoordinate.right *fullstate->low[2] );
    substate.high[2] = static_cast<short>((1-subcoordinate.right) *fullstate->high[0] + subcoordinate.right *fullstate->high[2]);
    substate.low[1] =  static_cast<short>((1-subcoordinate.top)   *fullstate->low[1]  + subcoordinate.top   *fullstate->low[3] );
    substate.high[1] = static_cast<short>((1-subcoordinate.top)   *fullstate->high[1] + subcoordinate.top   *fullstate->high[3]);
    substate.low[3] =  static_cast<short>((1-subcoordinate.bottom)*fullstate->low[1]  + subcoordinate.bottom*fullstate->low[3] );
    substate.high[3] = static_cast<short>((1-subcoordinate.bottom)*fullstate->high[1] + subcoordinate.bottom*fullstate->high[3]);
    return substate;  // by value
}

void PyramidQualityFunction::setup(int argnumpoints, int argwidth, int argheight, 
                               double* argxpos, double* argypos, double* argclst, 
                               void* argdata) {
    PyramidParameters* data = reinterpret_cast<PyramidParameters*>(argdata);
                               
    width = argwidth;
    height = argheight;

    cell_coordinates.clear();
    for (unsigned int l=1;l<=data->numlevels;l++) {
        for (unsigned int i=0;i<l;i++) {
            for (unsigned int j=0;j<l;j++) {
                Cell cur_cell;
                cur_cell.left = j/(float)l;
                cur_cell.top = i/(float)l;
                cur_cell.right = (j+1)/(float)l;
                cur_cell.bottom = (i+1)/(float)l;
                cell_coordinates.push_back(cur_cell);
            }
        }
    }
    unsigned int numcells = cell_coordinates.size();
    cell_weights.resize(numcells);
    for (unsigned int i=0; i<numcells; i++)
        cell_weights[i] = 1.;   // weighting comes later

    cell_quality.resize(numcells);
    for (unsigned int i=0; i<numcells; i++) {
        double* argweight = data->weightptr[i];
        cell_quality[i].setup(argnumpoints, argwidth, argheight, argxpos, argypos, argclst, argweight);
    }

    return;
}

void PyramidQualityFunction::cleanup() {
    return;
}

double PyramidQualityFunction::upper_bound(const sstate* state) const {
    double quality_bound=0.;
    for (unsigned int i=0; i<cell_quality.size(); i++) {
        sstate substate = rel_to_abs_coordinate(cell_coordinates[i], state);
        quality_bound += cell_weights[i] * cell_quality[i].upper_bound(&substate);
    }
    return quality_bound;
}

