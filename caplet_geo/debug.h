/*
CREATED : Jan 31, 2013
AUTHOR  : Yu-Chung Hsiao
EMAIL   : project.caplet@gmail.com

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


#ifndef DEBUG_H
#define DEBUG_H

#define DEBUG

//_________________________________________________________
//* geoloader
//#define DEBUG_PANEL_DISCRETIZATION
//#define DEBUG_CUT
//#define DEBUG_POLY2RECT



//_________________________________________________________
//* gdsgeometry
//#define CONDUCTOR_PRINT_ALL_RECTS
//- COMBINED_SHAPE
//  comment for more number of basis functions but better accuracy
//#define COMBINED_SHAPE

//_________________________________________________________
//* panelrenderer
#define WHITE_BACKGROUND

//_________________________________________________________
//* Parameters

//* Uncomment to use obsolete mergeProjection() algorithm
//* - Everything on the same support is merged when possible
//* - No projection distance info is used
//#define MERGE_PROJECTION_VER1_0

namespace caplet{
    const float DEFAULT_PROJECTION_MERGE_DISTANCE = 1e-7f;
    const float DEFAULT_PROJECTION_DISTANCE = 2e-6f;
    const float DEFAULT_COINCIDENTAL_MARGIN = 0.05;
};

#endif // DEBUG_H
