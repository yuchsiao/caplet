/*
Created : Jan 29, 2013
Modified: Feb 13, 2013
Author  : Yu-Chung Hsiao
Email   : yuchsiao@mit.edu
*/

/*
This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef CAPLET_CONST_H
#define CAPLET_CONST_H

namespace caplet{

enum DIR{
    X, Y, Z, FLAT
};
enum {
    MIN, MAX, LENGTH, CENTER
};

const int nDim = 3;
const int nBit = 4;
const float MAX_ASPECT_RATIO = 50;


typedef float (*shape_t)(float, float);

}

#endif // CAPLET_CONST_H
