/*
Created : Jan 29, 2013
Author  : Yu-Chung Hsiao
Email   : project.caplet@gmail.com
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


#ifndef CAPLET_PARAMETER_H
#define CAPLET_PARAMETER_H

//* Switches for using tabulated atan and log for each integral to speed up
//- Default: commented
//#define CAPLET_INIT_ATAN_LOG
//#define CAPLET_ATAN_LOG_INT_XY
//#define CAPLET_ATAN_LOG_INT_XYY
//#define CAPLET_ATAN_LOG_INT_XYZ
//#define CAPLET_ATAN_LOG_INT_XYXY
//#define CAPLET_ATAN_LOG_INT_XYYZ

//* Define to enable timer
#define CAPLET_TIMER
//* Define the number of executing extraction
//  in order to have averaged timing information
#define N_ITER 1

//* Use flat arch or side
//  Current arch and side functions are extracted from crossing wires
//  with 0.2um seperation and valid up to 0.4um
//  It may be sufficient for other separations too.
//  Using flat shapes for arches and sides may also be sufficient.
//- Default: uncommented (flat)
#define CAPLET_FLAT_ARCH
#define CAPLET_FLAT_SIDE

#define ROBUST_INTEGRAL_CHECK 

//* Openmp num of threads
#ifdef CAPLET_OPENMP
    #define CAPLET_OPENMP_NUM_THREADS 4
#endif

namespace caplet{

//* Gauss quad points and subdivision number setting
const int gauss_n = 2;

const int gauss_n_ZXZF = gauss_n;
const int gauss_n_ZXXF = gauss_n;
const int gauss_n_ZXYF = gauss_n;

const int gauss_n_ZXYX_1 = gauss_n;
const int gauss_n_ZXYX_2 = gauss_n;

const int gauss_n_ZXZY_1 = gauss_n;
const int gauss_n_ZXZY_2 = gauss_n;

const int gauss_n_ZXYZ_1 = gauss_n;
const int gauss_n_ZXYZ_2 = gauss_n;

const int gauss_n_ZXXY_1 = gauss_n;
const int gauss_n_ZXXY_2 = gauss_n;

const int quad_n = 2;

const int quad_n_ZXZX_1 = quad_n;
const int quad_n_ZXZX_2 = quad_n;

const int quad_n_ZXXZ_1 = quad_n;
const int quad_n_ZXXZ_2 = quad_n;

//* Switches for whether using fully analytical integrals or not
//  only work for all flat shapes
//- Default: all false;
const bool switch_analytical_integral = false;

const bool switch_analytical_intZXZF = switch_analytical_integral;
const bool switch_analytical_intZXXF = switch_analytical_integral;
const bool switch_analytical_intZXYF = switch_analytical_integral;

const bool switch_analytical_intZXZX = switch_analytical_integral;
const bool switch_analytical_intZXYX = switch_analytical_integral;
const bool switch_analytical_intZXZY = switch_analytical_integral;
const bool switch_analytical_intZXXZ = switch_analytical_integral;
const bool switch_analytical_intZXYZ = switch_analytical_integral;
const bool switch_analytical_intZXXY = switch_analytical_integral;


//* Define zero squared level used in analytical integrals
//  Consider the nano scale of VLSI interconnects and layout minimum grids.
//  Taking the length much lower than the grid size (nm) is a reasonable assumption
//- Default: 1e-12
const float zero  = 1e-12;
const float zero2 = zero*zero;


//* Approximation gaurd ring to increase guard ring for
//  approximating calColD by 1/r
const double approximationGuardRingForCalColD = 2.0;

}

#endif // CAPLET_PARAMETER_H
