/********************************************************
 *                                                      *
 *  Efficient Subwindow Search (ESS) implemented in C++ *
 *  generic class definintion for quality functions     *
 *                                                      *
 *   Copyright 2006-2008 Christoph Lampert              *
 *   Contact: <mail@christoph-lampert.org>              *
 *                                                      *
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

#ifndef _QUALITY_FUNCTION_H
#define _QUALITY_FUNCTION_H

class sstate;

class QualityFunction {
    public:

        // use the argdata field to pass in any additional data necessary
        virtual void setup(int argnumpoints, 
                                        int argwidth, 
                                        int argheight, 
                                        double* argxpos, 
                                        double* argypos, 
                                        double* argclst, 
                                        void* argdata) { return; };
                    
        virtual void cleanup() { return; };

        virtual double upper_bound(const sstate* state) const = 0;
};

#endif
