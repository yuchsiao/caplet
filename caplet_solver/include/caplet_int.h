/*
Created : Jul 31, 2010
Modified: Feb 14, 2013
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

#ifndef CAPLET_INT_H_
#define CAPLET_INT_H_

#include "caplet_debug.h"
#include "caplet_const.h"

#ifdef  DEBUG_SHAPE_BOUNDARY_CHECK
	#include <iostream>
#endif


namespace caplet{

//* FAST GALERKIN MODE

//* Flat-Flat integrals
//* - integral 1
float intZFZF(float* coord1[nDim][nBit], float* coord2[nDim][nBit]);
//* - integral 2
float intZFXF(float* coord1[nDim][nBit], float* coord2[nDim][nBit]);

//* Linear-Flat integrals
//* - integral 3
float intZXZF(
        float* coord1[nDim][nBit], float bz, float bsh, float (*shape)(float, float),
        float* coord2[nDim][nBit] );
//* - integral 4
float intZXXF(
        float* coord1[nDim][nBit], float bz, float bsh, float (*shape)(float, float),
        float* coord2[nDim][nBit] );
//* - integral 5
float intZXYF(
        float* coord1[nDim][nBit], float bz, float bsh, float (*shape)(float, float),
        float* coord2[nDim][nBit] );

//* Linear-Linear integrals
//* - integral 6
float intZXZX(
        float* coord1[nDim][nBit], float bz1, float bsh1, float (*shape1)(float, float),
        float* coord2[nDim][nBit], float bz2, float bsh2, float (*shape2)(float, float) );
//* - integral 7
float intZXYX(
        float* coord1[nDim][nBit], float bz1, float bsh1, float (*shape1)(float, float),
        float* coord2[nDim][nBit], float bz2, float bsh2, float (*shape2)(float, float) );
//* - integral 8
float intZXZY(
        float* coord1[nDim][nBit], float bz1, float bsh1, float (*shape1)(float, float),
        float* coord2[nDim][nBit], float bz2, float bsh2, float (*shape2)(float, float) );
//* - integral 9
float intZXXZ(
        float* coord1[nDim][nBit], float bz1, float bsh1, float (*shape1)(float, float),
        float* coord2[nDim][nBit], float bz2, float bsh2, float (*shape2)(float, float) );
//* - integral 10
float intZXYZ(
        float* coord1[nDim][nBit], float bz1, float bsh1, float (*shape1)(float, float),
        float* coord2[nDim][nBit], float bz2, float bsh2, float (*shape2)(float, float) );
//* - integral 11
float intZXXY(
        float* coord1[nDim][nBit], float bz1, float bsh1, float (*shape1)(float, float),
        float* coord2[nDim][nBit], float bz2, float bsh2, float (*shape2)(float, float) );



//***********************************************************
//*
//* Double collocation mode
//*
//*
double calColD(float* p1[nDim][nBit], float* p2[nDim][nBit]);



//**********************************
//*
//* Shape functions
//*
//*
inline float flat(float x, float w){
	return 1;
}


inline float arch(float x, float w){
    //* Constant arch shapes are good enough for most of the time
    #ifdef CAPLET_FLAT_ARCH
    return 0.5;
    #endif

    #ifdef DEBUG_SHAPE_BOUNDARY_CHECK
	if ( x>1.41e-6 || x<0 || w<0){
		std::cerr << "ERROR: input argument is out of the arch boundary" << std::endl;
	}
    #endif

    //* Extracted for arch length = 0.4um
    if ( w>10e-6 ){
		w = 10e-6;
	}
    return ( -1.1065e-06 + (241.98 + 5.1765e+06*w)*w + (350.58 + -2.3076e+08*x + -5.5044e+07*w)*x )/( -1.0785e-05 + ( 404.46 + ( 2.6513e+06 + 2.7367e+11*w )*w )*w + ( 456.37 + (1.8301e+08 + -3.3833e+14*x)*x + (-1.0181e+07 + 5.4154e+12*w + 1.0317e+15*x )*w )*x  );
}


inline float side(float x, float w){
    #ifdef CAPLET_FLAT_SIDE
    return 0.5;
    #endif

    #ifdef DEBUG_SHAPE_BOUNDARY_CHECK
	if ( x>1.251e-6 || x<0 || w<0){
		std::cerr << "ERROR: input argument is out of the arch boundary" << std::endl;
	}
    #endif

    //* Extracted for arch length = 0.4um
    if ( w>10e-6 ){
		w = 10e-6;
	}
	return ( 0.00012065 + ( 237.09 + -8.6825e+06*w )*w + ( 1580.3 + -1.7957e+09*x + -2.3054e+07*w )*x )/( 0.00029255 + ( 398.95 + -1.6467e+07*w )*w + ( 2789.4 + 1.447e+09*x + 1.0577e+08*w )*x );
}


}

#endif /* CAPLET_INT_H_ */
