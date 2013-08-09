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

#include "caplet_int.h"


#include "caplet_parameter.h"
#include "caplet_elem.h"
#include "caplet_blas.h"
#include "caplet_widgets.h"
#include "caplet_gauss.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>


namespace caplet{

using std::abs;
using std::sqrt;


//* Integrals
//- subroutine declaration
float  int_xyxy(float* p1[nDim][nBit], float* p2[nDim][nBit]);

float  int_xyyz(float a, float b, float ly, float lz, float x, float y, float z);
float  int_xyyz(float* p1[nDim][nBit], float* p2[nDim][nBit]);
float  int_xyxy(float a, float b, float lx, float ly, float x, float y, float z);
float  int_xyzx(float* p1[nDim][nBit], float* p2[nDim][nBit]);

float  int_xyy(float a, float b, float ly, float x, float y, float z);
float  int_xyz(float a, float b, float lz, float x, float y, float z);

float  int_xy(float a, float b, float x, float y, float z, float area);



//* Three internal guass quad
//* 1. x-dir quad over int_xy
inline float gauss_int_xy_x(int gauss_n, float a, float b, float x1, float x2, float y, float z, float (*shape)(float, float), float w, float bz){
	int gauss_n2 = (gauss_n+1)/2;
	float pm = (x2+x1)/2;
	float pr = (x2-x1)/2;

	float  val=0;
	float  p0 = (bz>0)? (x1) : (x2);
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * int_xy(a, b, abs(pm), y, z, a*b) * shape(abs(pm-p0), w);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+int_xy(a, b, abs(pm-dp), y, z, a*b) * shape(abs(pm-dp-p0), w)
				+int_xy(a, b, abs(pm+dp), y, z, a*b) * shape(abs(pm+dp-p0), w)
				);
	}
	return val * pr;
}

//* 2. y-dir quad over int_xy
inline float gauss_int_xy_y(int gauss_n, float a, float b, float x, float y1, float y2, float z, float (*shape)(float, float), float w, float bz){
	int gauss_n2 = (gauss_n+1)/2;
	float pm = (y2+y1)/2;
	float pr = (y2-y1)/2;

	float  val=0;
	float  p0 = (bz>0)? (y1) : (y2);
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * int_xy(a, b, x, abs(pm), z, a*b) * shape(abs(pm-p0), w);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+int_xy(a, b, x, abs(pm-dp), z, a*b) * shape(abs(pm-dp-p0), w)
				+int_xy(a, b, x, abs(pm+dp), z, a*b) * shape(abs(pm+dp-p0), w)
				);
	}
	return val * pr;
}

//* 3. z-dir quad over int_xy
inline float gauss_int_xy_z(int gauss_n, float a, float b, float x, float y, float z1, float z2, float (*shape)(float , float ), float w, float bz){
	int gauss_n2 = (gauss_n+1)/2;
	float pm = (z2+z1)/2;
	float pr = (z2-z1)/2;

	float  val=0;
	float  p0 = (bz>0)? (z1) : (z2);
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * int_xy(a,b,x,y, abs(pm),a*b) * shape(abs(pm-p0), w);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+int_xy(a,b,x,y, abs(pm+dp),a*b) * shape(abs(pm+dp-p0), w)
				+int_xy(a,b,x,y, abs(pm-dp),a*b) * shape(abs(pm-dp-p0), w)
				);
	}
	return val * pr;
}



//* Flat-Flat integrals
//  - integral 1
float intZFZF(float* coord1[3][4], float* coord2[3][4]){
	return int_xyxy( coord1, coord2 );
}


//* - integral 2
float intZFXF(float* coord1[3][4], float* coord2[3][4]){
	return int_xyyz( coord1, coord2 );
}


//* Linear-Flat integrals
//* - integral 3
float intZXZF(
		float* coord1[3][4], float bz, float bsh, float (*shape)(float, float),
		float* coord2[3][4]
){
    //****
    if (switch_analytical_intZXZF == true){
        return int_xyxy( coord1, coord2 );
    }
    //****

	float a  = (*coord2[X])[LENGTH];
	float b  = (*coord2[Y])[LENGTH];
	float ly = (*coord1[Y])[LENGTH];
	float w  = ly;

	if ( b<ly ){
		float temp = b; b=ly; ly=temp;
	}
	float y = abs( (*coord2[Y])[CENTER] - (*coord1[Y])[CENTER] );
	float z = abs( (*coord2[Z])[CENTER] - (*coord1[Z])[CENTER] );

	float x1 = (*coord1[X])[0] - (*coord2[X])[CENTER];
	float x2 = (*coord1[X])[1] - (*coord2[X])[CENTER];


	//BEGIN____________________________gauss quad_______________________________
	const int gauss_n = gauss_n_ZXZF;

	const int gauss_n2 = (gauss_n+1)/2;
	const float pm = (x2+x1)/2;
	const float pr = (x2-x1)/2;
	const float  p0 = (bz>0)? (x1) : (x2);

	float  val=0;
	int init_i = 0;
	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * int_xyy(a,b,ly,abs(pm),y,z) * shape(abs(pm-p0), w);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
        float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+int_xyy(a,b,ly,abs(pm+dp),y,z) * shape(abs(pm+dp-p0), w)
				+int_xyy(a,b,ly,abs(pm-dp),y,z) * shape(abs(pm-dp-p0), w)
				);
	}
	return val * pr;
	//_END_____________________________gauss quad_______________________________
}


//* - integral 4
float intZXXF(
		float* coord1[3][4], float bz, float bsh, float (*shape)(float, float),
		float* coord2[3][4]
){
    //****
    if (switch_analytical_intZXXF == true){
        return int_xyyz( coord1, coord2 );
    }
    //****

	float a = (*coord2[Z])[LENGTH];
	float b = (*coord2[Y])[LENGTH];
	float ly = (*coord1[Y])[LENGTH];
	float w  = ly;

	if ( b<ly ){
		float temp = b; b = ly; ly = temp;
	}
	float x = abs( (*coord2[Z])[CENTER] - (*coord1[Z])[CENTER] );
	float y = abs( (*coord2[Y])[CENTER] - (*coord1[Y])[CENTER] );

	float x1 = (*coord1[X])[0] - (*coord2[X])[CENTER];
	float x2 = (*coord1[X])[1] - (*coord2[X])[CENTER];


	//BEGIN____________________________gauss quad_______________________________
	const int gauss_n = gauss_n_ZXXF;

	const int gauss_n2 = (gauss_n+1)/2;
	const float pm = (x2+x1)/2;
	const float pr = (x2-x1)/2;
	const float  p0 = (bz>0)? (x1) : (x2);

	float  val=0;
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * int_xyy(a,b,ly,x,y,abs(pm)) * shape(abs(pm-p0), w);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+int_xyy(a,b,ly,x,y, abs(pm+dp)) * shape(abs(pm+dp-p0), w)
				+int_xyy(a,b,ly,x,y, abs(pm-dp)) * shape(abs(pm-dp-p0), w)
				);
	}
	return val * pr;
	//_END_____________________________gauss quad_______________________________
}


//* - integral 5
float intZXYF(
		float* coord1[3][4], float bz, float bsh, float (*shape)(float, float),
		float* coord2[3][4]
){
    //****
    if (switch_analytical_intZXYF == true){
        return int_xyzx( coord1, coord2 );
    }
    //****

	float a = (*coord2[Z])[LENGTH];
	float b = (*coord2[X])[LENGTH];
	float lz= (*coord1[Y])[LENGTH]; // length in rotated z

	float w = lz;

	float x = abs((*coord2[Z])[CENTER] - (*coord1[Z])[CENTER]);
	float z = abs((*coord2[Y])[CENTER] - (*coord1[Y])[CENTER]);

	float y1 = (*coord1[X])[0] - (*coord2[X])[CENTER];
	float y2 = (*coord1[X])[1] - (*coord2[X])[CENTER];

	//BEGIN____________________________gauss quad_______________________________
	const int gauss_n = gauss_n_ZXYF;

	const int gauss_n2 = (gauss_n+1)/2;
	const float pm = (y2+y1)/2;
	const float pr = (y2-y1)/2;
	const float p0 = (bz>0)? (y1) : (y2);

	float  val=0;
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * int_xyz(a,b,lz,x,abs(pm),z) * shape(abs(pm-p0), w);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+int_xyz(a,b,lz,x,abs(pm+dp),z) * shape(abs(pm+dp-p0), w)
				+int_xyz(a,b,lz,x,abs(pm-dp),z) * shape(abs(pm-dp-p0), w)
				);
	}
	return val * pr;
	//_END_____________________________gauss quad_______________________________
}


//* Linear-Linear integrals
//* - integral 6
float intZXZX(
		float* coord1[3][4], float bz1, float bsh1, float (*shape1)(float, float),
		float* coord2[3][4], float bz2, float bsh2, float (*shape2)(float, float)
){
    //****
    if (switch_analytical_intZXZX == true){
        return int_xyxy( coord1, coord2 );
    }
    //****

	float w1 = (*coord1[Y])[LENGTH];
	float w2 = (*coord2[Y])[LENGTH];

	int quad_n1 = quad_n_ZXZX_1;
	int quad_n2 = quad_n_ZXZX_2;

	float a1 = (*coord1[X])[LENGTH] / quad_n1;
	float a2 = (*coord2[X])[LENGTH] / quad_n2;

	float val = 0;

	float z = abs( (*coord1[Z])[CENTER] - (*coord2[Z])[CENTER] );
	float y = abs( (*coord1[Y])[CENTER] - (*coord2[Y])[CENTER] );

	float x10 = (bz1>0)? ( (*coord1[X])[0] ) : ( (*coord1[X])[1] );
	float x20 = (bz2>0)? ( (*coord2[X])[0] ) : ( (*coord2[X])[1] );

	float x2 = (*coord2[X])[0] - a2/2;
	for ( int i=0; i<quad_n1; i++ ){
		x2 += a2;
		float temp = 0;
		float x1 = (*coord1[X])[0] - a1/2;
		for ( int j=0; j<quad_n2; j++ ){
			x1 += a1;
			float x = abs( x1-x2 );
			temp += int_xyxy(a1,w1,a2,w2,x,y,z)*shape1( abs(x1-x10), w1 );
		}
		val += temp*shape2( abs(x2-x20), w2 );
	}
	return val;
}


//* - integral 7
float intZXYX(
		float* coord1[3][4], float bz1, float bsh1, float (*shape1)(float, float),
		float* coord2[3][4], float bz2, float bsh2, float (*shape2)(float, float)
){
    //****
    if (switch_analytical_intZXYX == true){
        return int_xyzx( coord1, coord2 );
    }
    //****

	float a = (*coord1[Y])[LENGTH];
	float b = (*coord2[Z])[LENGTH];
	float w1 = a;
	float w2 = b;

	float x = abs((*coord2[Y])[CENTER]-(*coord1[Y])[CENTER]);
	float y = abs((*coord2[Z])[CENTER]-(*coord1[Z])[CENTER]);


	// output integration limits
	float zp1 = (*coord1[X])[0];
	float zp2 = (*coord1[X])[1];

	// inner integration limits
	float z1 = (*coord2[X])[0];
	float z2 = (*coord2[X])[1];


	//BEGIN____________________________gauss quad_______________________________
	const int gauss_n_inner = gauss_n_ZXYX_1;
	const int gauss_n = gauss_n_ZXYX_2;

	const int gauss_n2 = (gauss_n+1)/2;
	const float pm = (zp2+zp1)/2;
	const float pr = (zp2-zp1)/2;
	const float  p0 = (bz1>0)? (zp1) : (zp2);

	float  val=0;
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * gauss_int_xy_z(gauss_n_inner, a,b,x,y, z1-pm, z2-pm, shape2, w2, bz2) * shape1(abs(pm-p0), w1);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+gauss_int_xy_z(gauss_n_inner, a,b,x,y, z1-pm-dp, z2-pm-dp, shape2, w2, bz2) * shape1(abs(pm+dp-p0), w1)
				+gauss_int_xy_z(gauss_n_inner, a,b,x,y, z1-pm+dp, z2-pm+dp, shape2, w2, bz2) * shape1(abs(pm-dp-p0), w1)
				);
	}
	return val * pr;
	//_END_____________________________gauss quad_______________________________
}


//* - integral 8
float intZXZY(
		float* coord1[3][4], float bz1, float bsh1, float (*shape1)(float, float),
		float* coord2[3][4], float bz2, float bsh2, float (*shape2)(float, float)
){
    //****
    if (switch_analytical_intZXZY == true){
        return int_xyxy( coord1, coord2 );
    }
    //****

	float a = (*coord2[X])[LENGTH];
	float b = (*coord1[Y])[LENGTH];
	float w1 = b;
	float w2 = a;

	float z = abs((*coord2[Z])[CENTER] - (*coord1[Z])[CENTER]);

	// outer integration limits
	float x1 = (*coord1[X])[0] - (*coord2[X])[CENTER];
	float x2 = (*coord1[X])[1] - (*coord2[X])[CENTER];


	// inner integration limits
	float y1 = (*coord2[Y])[0] - (*coord1[Y])[CENTER];
	float y2 = (*coord2[Y])[1] - (*coord1[Y])[CENTER];


	//BEGIN____________________________gauss quad_______________________________
	const int gauss_n_inner = gauss_n_ZXZY_1;
	const int gauss_n = gauss_n_ZXZY_2;

	const int gauss_n2 = (gauss_n+1)/2;
	const float pm = (x2+x1)/2;
	const float pr = (x2-x1)/2;
	const float  p0 = (bz1>0)? (x1) : (x2);

	float  val=0;
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * gauss_int_xy_y(gauss_n_inner, a,b, abs(pm), y1, y2, z, shape2, w2, bz2) * shape1(abs(pm-p0), w1);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+gauss_int_xy_y(gauss_n_inner, a,b, abs(pm-dp), y1, y2, z, shape2, w2, bz2) * shape1(abs(pm-dp-p0), w1)
				+gauss_int_xy_y(gauss_n_inner, a,b, abs(pm+dp), y1, y2, z, shape2, w2, bz2) * shape1(abs(pm+dp-p0), w1)
				);
	}
	return val * pr;
	//_END_____________________________gauss quad_______________________________
}


//* - integral 9
float intZXXZ(
		float* coord1[3][4], float bz1, float bsh1, float (*shape1)(float, float),
		float* coord2[3][4], float bz2, float bsh2, float (*shape2)(float, float)
){
    //****
    if (switch_analytical_intZXXZ == true){
        return int_xyyz( coord1, coord2 );
    }
    //****

	float w1 = (*coord1[Y])[LENGTH];
	float w2 = (*coord2[Y])[LENGTH];

	int quad_n1 = quad_n_ZXXZ_1;
	int quad_n2 = quad_n_ZXXZ_2;

	float a1 = (*coord1[X])[LENGTH] / quad_n1;
	float lz = (*coord2[Z])[LENGTH] / quad_n2;

	float y = abs( (*coord1[Y])[CENTER] - (*coord2[Y])[CENTER] );
	float x2 = (*coord2[X])[CENTER];
	float z1 = (*coord1[Z])[CENTER];

	float x0 = (bz1>0)? ( (*coord1[X])[0] ) : ( (*coord1[X])[1] );
	float z0 = (bz2>0)? ( (*coord2[Z])[0] ) : ( (*coord2[Z])[1] );

	float val = 0;
	float z2 = (*coord2[Z])[0] - lz/2;
	for ( int i=0; i<quad_n2; i++ ){
		z2 += lz;
		float temp = 0;
		float x1 = (*coord1[X])[0] - a1/2;
		for ( int j=0; j<quad_n1; j++ ){
			x1 += a1;
			float x = abs( x1 - x2 );
			float z = abs( z2 - z1 );
			temp += int_xyyz(a1,w1,w2,lz,x,y,z)*shape1( abs(x1-x0), w1 );
		}
		val += temp*shape2( abs(z2-z0), w2 );
	}
	return val;
}


//* - integral 10
float intZXYZ(
		float* coord1[3][4], float bz1, float bsh1, float (*shape1)(float, float),
		float* coord2[3][4], float bz2, float bsh2, float (*shape2)(float, float)
){
    //****
    if (switch_analytical_intZXYZ == true){
        return int_xyzx( coord1, coord2 );
    }
    //****

	float a = (*coord2[X])[LENGTH];
	float b = (*coord1[Y])[LENGTH];
	float w1 = b;
	float w2 = a;

	float y = abs((*coord2[Y])[CENTER] - (*coord1[Y])[CENTER]);

	// outer integration limits
	float x1 = (*coord1[X])[0] - (*coord2[X])[CENTER];
	float x2 = (*coord1[X])[1] - (*coord2[X])[CENTER];

	// inner integration limits
	float z1 = (*coord2[Z])[0] - (*coord1[Z])[CENTER];
	float z2 = (*coord2[Z])[1] - (*coord1[Z])[CENTER];

	//BEGIN____________________________gauss quad_______________________________
	const int gauss_n_inner = gauss_n_ZXYZ_1;
	const int gauss_n = gauss_n_ZXYZ_2;

	const int gauss_n2 = (gauss_n+1)/2;
	const float pm = (x2+x1)/2;
	const float pr = (x2-x1)/2;
	const float  p0 = (bz1>0)? (x1) : (x2);

	float  val=0;
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * gauss_int_xy_z(gauss_n_inner, a,b, abs(pm) ,y, z1, z2, shape2, w2, bz2) * shape1(abs(pm-p0), w1);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+gauss_int_xy_z(gauss_n_inner, a,b, abs(pm+dp) ,y, z1, z2, shape2, w2, bz2) * shape1(abs(pm+dp-p0), w1)
				+gauss_int_xy_z(gauss_n_inner, a,b, abs(pm-dp) ,y, z1, z2, shape2, w2, bz2) * shape1(abs(pm-dp-p0), w1)
				);
	}
	return val * pr;
	//_END_____________________________gauss quad_______________________________
}


//* - integral 11
float intZXXY(
		float* coord1[3][4], float bz1, float bsh1, float (*shape1)(float, float),
		float* coord2[3][4], float bz2, float bsh2, float (*shape2)(float, float)
){
    //****
    if (switch_analytical_intZXXY == true){
        return int_xyyz( coord1, coord2 );
    }
    //****

	float a = (*coord1[Y])[LENGTH];
	float b = (*coord2[Z])[LENGTH];
	float w1 = a;
	float w2 = b;

	float y = abs((*coord2[Z])[CENTER] - (*coord1[Z])[CENTER]);

	// outer intergration limits
	float z1 = (*coord1[X])[0] - (*coord2[X])[CENTER];
	float z2 = (*coord1[X])[1] - (*coord2[X])[CENTER];

	// inner integration limits
	float x1 = (*coord2[Y])[0]-(*coord1[Y])[CENTER];
	float x2 = (*coord2[Y])[1]-(*coord1[Y])[CENTER];

	//BEGIN____________________________gauss quad_______________________________
	const int gauss_n_inner = gauss_n_ZXXY_1;
	const int gauss_n = gauss_n_ZXXY_2;

	const int gauss_n2 = (gauss_n+1)/2;
	const float pm = (z2+z1)/2;
	const float pr = (z2-z1)/2;
	const float  p0 = (bz1>0)? (z1) : (z2);

	float  val=0;
	int init_i = 0;

	if ( (gauss_n%2)==1 ){ // odd n
		init_i = 1;
		val += (*gauss::w[gauss_n])[0] * gauss_int_xy_x(gauss_n_inner, a,b, x1, x2, y, abs(pm), shape2, w2, bz2) * shape1(abs(pm-p0), w1);
	}
	for ( int i = init_i; i < gauss_n2; i++ ) {
		float dp = pr * (*gauss::p[gauss_n])[i];
		val += (*gauss::w[gauss_n])[i] * (
				+gauss_int_xy_x(gauss_n_inner, a,b, x1, x2, y, abs(pm-dp), shape2, w2, bz2) * shape1(abs(pm-dp-p0), w1)
				+gauss_int_xy_x(gauss_n_inner, a,b, x1, x2, y, abs(pm+dp), shape2, w2, bz2) * shape1(abs(pm+dp-p0), w1)
				);
	}
	return val * pr;
	//_END_____________________________gauss quad_______________________________
}





//***************************************
//*
//* Internal integration
//*
//*

float  int_xy(float a, float b, float x, float y, float z, float area){
    //**** ASSERT
    assert( a>0 );
	assert( b>0 );
	assert( x>=0 );
	assert( y>=0 );
	assert( z>=0 );
	assert( area>0 );
    //****

    //**** Select atan log verion
    #ifdef CAPLET_ATAN_LOG_INT_XY
    using caplet::atan;
    using caplet::log;
    #else
    using std::atan;
    using std::log;
    #endif
    //****

    //* Compute aspect ratio
	float ba = b/a;

	if ( ba < 1 ){ // swap x and y
		ba = 1/ba;
		float temp;
		temp = b; b = a; a = temp;
		temp = y; y = x; x = temp;
	}

    //* Approximation when the evaluation point is far from the origin
	if ( ba < 8 ){
		if ( x > 4.545f*b
				 ||
			 y > ( -0.036938f+0.044672f*ba )/( -0.0051053f+0.0068095f*ba ) *b
			 	 ||
			 z > ( -0.011618f+0.029056f*ba )/( -0.0037776f+0.0064776f*ba ) *b
		){
            return area/sqrt(x*x+y*y+z*z);
		}
	}else{
		if ( x > 4.545f*b || y > 6.545*b || z > 4.636*b ){
            return area/sqrt(x*x+y*y+z*z);
		}
	}

    //* Analytical integral
    float x1 = x-a/2;
    float x2 = x+a/2;
    float y1 = y-b/2;
    float y2 = y+b/2;

    float x1_2 = x1*x1;
    float x2_2 = x2*x2;
    float y1_2 = y1*y1;
    float y2_2 = y2*y2;

    float z_2  = z*z;

    float r11 = sqrt(x1_2 + y1_2 + z_2);
    float r12 = sqrt(x1_2 + y2_2 + z_2);
    float r21 = sqrt(x2_2 + y1_2 + z_2);
    float r22 = sqrt(x2_2 + y2_2 + z_2);

    #ifdef ROBUST_INTEGRAL_CHECK
    float p =  (( x2_2>zero2 && abs(y2+r22)>zero && abs(y1+r21)>zero )? x2*log( (y2+r22)/(y1+r21) ) :0)
              +(( x1_2>zero2 && abs(y1+r11)>zero && abs(y2+r12)>zero )? x1*log( (y1+r11)/(y2+r12) ) :0)
              +(( y2_2>zero2 && abs(x2+r22)>zero && abs(x1+r12)>zero )? y2*log( (x2+r22)/(x1+r12) ) :0)
              +(( y1_2>zero2 && abs(x1+r11)>zero && abs(x2+r21)>zero )? y1*log( (x1+r11)/(x2+r21) ) :0)
              +((  z_2>zero2 )?  z*( atan(x2*y1/r21/z) + atan(x1*y2/r12/z)
                                    -atan(x2*y2/r22/z) - atan(x1*y1/r11/z) ):0);
    #else
    float p =  ((x2_2>zero2)? x2*log( (y2+r22)/(y1+r21) ):0)
              +((x1_2>zero2)? x1*log( (y1+r11)/(y2+r12) ):0)
              +((y2_2>zero2)? y2*log( (x2+r22)/(x1+r12) ):0)
              +((y1_2>zero2)? y1*log( (x1+r11)/(x2+r21) ):0)
              +(( z_2>zero2)?  z*( atan(x2*y1/r21/z) + atan(x1*y2/r12/z)
                                  -atan(x2*y2/r22/z) - atan(x1*y1/r11/z) ):0);
    #endif
    return p;
}


float int_xyxy(float a, float b, float lx, float ly, float x, float y, float z){
    float xp[nBit];
    float yp[nBit];
    float zp[nBit];

    float xx[nBit];
    float yy[nBit];
    float zz[nBit];

    xp[MIN] = -a/2;
    xp[MAX] =  a/2;
    yp[MIN] = -b/2;
    yp[MAX] =  b/2;
    zp[MIN] =  0;
    zp[MAX] =  0;

	xp[CENTER] = 0;
	yp[CENTER] = 0;
	zp[CENTER] = 0;

	xp[LENGTH] = a;
	yp[LENGTH] = b;
	zp[LENGTH] = 0;


    xx[MIN] = -lx/2 +x;
    xx[MAX] =  lx/2 +x;
    yy[MIN] = -ly/2 +y;
    yy[MAX] =  ly/2 +y;
    zz[MIN] =  z;
    zz[MAX] =  z;

	xx[CENTER] = x;
	yy[CENTER] = y;
	zz[CENTER] = z;

	xx[LENGTH] = lx;
	yy[LENGTH] = ly;
	zz[LENGTH] = 0;

    float* coord1[nDim][nBit];
    float* coord2[nDim][nBit];

    *coord1[X] = xp;
    *coord1[Y] = yp;
    *coord1[Z] = zp;

    *coord2[X] = xx;
    *coord2[Y] = yy;
    *coord2[Z] = zz;

    return int_xyxy(coord1, coord2);
}


float int_xyxy(float* p1[3][4], float* p2[3][4]){
    //* note: the precision needs to be 'double' due to accuracy

    #ifdef CAPLET_ATAN_LOG_INT_XYXY
    using caplet::atan;
    using caplet::log;
    #else
    using std::atan;
    using std::log;
    #endif

	const float* xp;
	const float* yp;

	const float* xx;
	const float* yy;

    //* Assign longer sides in x- and y-dir to panel'(xp, yp)
    //  so that panel(xx, yy) is smaller on both sides
	if ( (*p1[X])[LENGTH] > (*p2[X])[LENGTH] ){
		xp = *p1[X];
		xx = *p2[X];
	}else{
		xp = *p2[X];
		xx = *p1[X];
	}
	if ( (*p1[Y])[LENGTH] > (*p2[Y])[LENGTH] ){
		yp = *p1[Y];
		yy = *p2[Y];
	}else{
		yp = *p2[Y];
		yy = *p1[Y];
	}

    //* Now xp.length > xx.length and yp.length > yy.length
    //  Next, determine which one is the characteristic length
	float b, lx, ly, a, x, y;
	if ( yp[LENGTH] > xp[LENGTH] ){
		a  = xp[LENGTH];
		b  = yp[LENGTH];
		lx = xx[LENGTH];
		ly = yy[LENGTH];
        //* calculate xcb,ycb,z
		x = abs(xp[CENTER]-xx[CENTER]);
		y = abs(yp[CENTER]-yy[CENTER]);
	}else{
		a  = yp[LENGTH];
		b  = xp[LENGTH];
		lx = yy[LENGTH];
		ly = xx[LENGTH];
        //* calculate xcb,ycb,z
		x = abs(yp[CENTER]-yy[CENTER]);
		y = abs(xp[CENTER]-xx[CENTER]);
	}
	float ba  = b/a;
    float lyb = ly/b;	// [0.001, 1]

	float z  = abs( (*p1[Z])[CENTER] - (*p2[Z])[CENTER] );


    //* Approximation when the separation between two panels is far enough
    //  Different approximation formula for different aspect ratios of panel'
    if ( ba>=6 ){
        if ( x > (( lyb<0.3 )? 0.3*b : ( -0.84064 + ( 4.2879 - 1.2213 *lyb )* lyb )*b )
                ||
             y > (( lyb<0.3 )? b     : ( 0.25582 + 3.0718*lyb )*b )
                ||
             z > ( 0.5041 -0.0053355*ba +2.0298*lyb  )*b
        ){
            return int_xy( a, b, x, y, z, a*b )*lx*ly;
        }
    }else if( ba>=4 ){
        if (( x > (( lyb<0.3 )? 0.3*b : ( -0.84064 + ( 4.2879 - 1.2213 *lyb )* lyb )*b ))
                ||
            ( y > (( lyb<0.3 )? b     : ( 0.25582 + 3.0718*lyb )*b ))
                ||
            ( z > ( ( 0.0028979-0.00023193*ba+0.0035999*lyb)/(-2.0903e-05+ ( 0.0013204 -9.0076e-05*ba )*ba + ( 0.002485 -0.0010111*lyb -0.00070171*ba )*lyb ) )*b)
        ){

            return int_xy( a, b, x, y, z, a*b )*lx*ly;
        }
    }else{ // ba<4
        if ( x > ( 3 )*b
                ||
             y > (  3.3157+ (-1.2913 + 0.17667*ba )*ba + (-3.2682 + 3.0242 * lyb + 0.79055*ba )*lyb   ) *b
                ||
             z > ( ( 0.0028979-0.00023193*ba+0.0035999*lyb)/(-2.0903e-05+ ( 0.0013204 -9.0076e-05*ba )*ba + ( 0.002485 -0.0010111*lyb -0.00070171*ba )*lyb ) )*b
         ){
            return int_xy( a, b, x, y, z, a*b )*lx*ly;
         }
    }

    //* Analytical integral
    //    double value 	 = 0.0f;
    //    const  double z2 = z*z;
    //    double u, u2, v, v2, r, f11;

    //    for ( int m=0; m<2; m++ ){
    //        for ( int n=0; n<2; n++ ){
    //            u = xp[m]-xx[n];
    //            u2= u*u;

    //            for ( int p=0; p<2; p++ ){
    //                for ( int q=0; q<2; q++ ){
    //                    v = yp[p]-yy[q];
    //                    v2= v*v;
    //                    r = sqrt(u2+v2+z2);

    //                    f11 = 0;
    //                    if ( z2 > zero2 ){
    //                        f11 +=
    //                        -12*u*v*z*atan(u*v/r/z)
    //                        +4*r*z2
    //                        +6*v*z2*log(r-v)
    //                        +6*u*z2*log(r-u);
    //                    }
    //                    if ( u2 > zero2 ){
    //                        f11 +=
    //                                -2*r*u2
    //                                -3*u2*v
    //                                +6*u2*v*log(r+v);
    //                    }
    //                    if ( v2 > zero2 ){
    //                        f11 +=
    //                                -2*r*v2
    //                                -3*u*v2
    //                                +6*u*v2*log(r+u);
    //                    }

    //                    value += (1-2*((m+n+p+q)%2))*f11;
    //                }
    //            }
    //        }
    //    }
    //    return value/12;

    //* Expanded analytical integral to speed up execution
    const double z2 	= z*z;

    const double u00		= xp[0] - xx[0];		const double u002 = u00 * u00;
    const double u01		= xp[0] - xx[1];		const double u012 = u01 * u01;
    const double u10		= xp[1] - xx[0];		const double u102 = u10 * u10;
    const double u11		= xp[1] - xx[1];		const double u112 = u11 * u11;

    const double v00		= yp[0] - yy[0];		const double v002 = v00 * v00;
    const double v01		= yp[0] - yy[1];		const double v012 = v01 * v01;
    const double v10		= yp[1] - yy[0];		const double v102 = v10 * v10;
    const double v11		= yp[1] - yy[1];		const double v112 = v11 * v11;

    const double r0000	= sqrt( u002 + v002 + z2);
    const double r0001	= sqrt( u002 + v012 + z2);
    const double r0010	= sqrt( u002 + v102 + z2);
    const double r0011	= sqrt( u002 + v112 + z2);
    const double r0100	= sqrt( u012 + v002 + z2);
    const double r0101	= sqrt( u012 + v012 + z2);
    const double r0110	= sqrt( u012 + v102 + z2);
    const double r0111	= sqrt( u012 + v112 + z2);
    const double r1000	= sqrt( u102 + v002 + z2);
    const double r1001	= sqrt( u102 + v012 + z2);
    const double r1010	= sqrt( u102 + v102 + z2);
    const double r1011	= sqrt( u102 + v112 + z2);
    const double r1100	= sqrt( u112 + v002 + z2);
    const double r1101	= sqrt( u112 + v012 + z2);
    const double r1110	= sqrt( u112 + v102 + z2);
    const double r1111	= sqrt( u112 + v112 + z2);

    const double u00v00 = u00 * v00;
    const double u00v01 = u00 * v01;
    const double u00v10 = u00 * v10;
    const double u00v11 = u00 * v11;

    const double u01v00 = u01 * v00;
    const double u01v01 = u01 * v01;
    const double u01v10 = u01 * v10;
    const double u01v11 = u01 * v11;

    const double u10v00 = u10 * v00;
    const double u10v01 = u10 * v01;
    const double u10v10 = u10 * v10;
    const double u10v11 = u10 * v11;

    const double u11v00 = u11 * v00;
    const double u11v01 = u11 * v01;
    const double u11v10 = u11 * v10;
    const double u11v11 = u11 * v11;


    double value =
            -2*(+r0000 *( u002 + v002 )
                -r0001 *( u002 + v012 )
                -r0010 *( u002 + v102 )
                +r0011 *( u002 + v112 )

                -r0100 *( u012 + v002 )
                +r0101 *( u012 + v012 )
                +r0110 *( u012 + v102 )
                -r0111 *( u012 + v112 )

                -r1000 *( u102 + v002 )
                +r1001 *( u102 + v012 )
                +r1010 *( u102 + v102 )
                -r1011 *( u102 + v112 )

                +r1100 *( u112 + v002 )
                -r1101 *( u112 + v012 )
                -r1110 *( u112 + v102 )
                +r1111 *( u112 + v112 ) );
    double temp = 0;

    #ifdef ROBUST_INTEGRAL_CHECK
    if (u002>zero2 &&
        abs(r0000+v00)>zero &&
        abs(r0001+v01)>zero &&
        abs(r0010+v10)>zero &&
        abs(r0011+v11)>zero    ){
    #else
    if (u002>zero2){
    #endif
        temp +=  u00*(	+u00v00*(-1 + 2*log(r0000+v00))
                        -u00v01*(-1 + 2*log(r0001+v01))
                        -u00v10*(-1 + 2*log(r0010+v10))
                        +u00v11*(-1 + 2*log(r0011+v11))	);
    }

    #ifdef ROBUST_INTEGRAL_CHECK
    if (u012>zero2 &&
        abs(r0100+v00)>zero &&
        abs(r0101+v01)>zero &&
        abs(r0110+v10)>zero &&
        abs(r0111+v11)>zero    ){
    #else
    if (u012>zero2){
    #endif
        temp += -u01*( +u01v00*(-1 + 2*log(r0100+v00))
                       -u01v01*(-1 + 2*log(r0101+v01))
                       -u01v10*(-1 + 2*log(r0110+v10))
                       +u01v11*(-1 + 2*log(r0111+v11)) );

    }

    #ifdef ROBUST_INTEGRAL_CHECK
    if (u102>zero2 &&
        abs(r1000+v00)>zero &&
        abs(r1001+v01)>zero &&
        abs(r1010+v10)>zero &&
        abs(r1011+v11)>zero    ){
    #else
    if (u102>zero2){
    #endif
        temp += -u10*( +u10v00*(-1 + 2*log(r1000+v00))
                       -u10v01*(-1 + 2*log(r1001+v01))
                       -u10v10*(-1 + 2*log(r1010+v10))
                       +u10v11*(-1 + 2*log(r1011+v11)) );
    }

    #ifdef ROBUST_INTEGRAL_CHECK
    if (u112>zero2 &&
        abs(r1100+v00)>zero &&
        abs(r1101+v01)>zero &&
        abs(r1110+v10)>zero &&
        abs(r1111+v11)>zero    ){
    #else
    if (u112>zero2){
    #endif
        temp += +u11*( +u11v00*(-1 + 2*log(r1100+v00))
                       -u11v01*(-1 + 2*log(r1101+v01))
                       -u11v10*(-1 + 2*log(r1110+v10))
                       +u11v11*(-1 + 2*log(r1111+v11)) );
    }

    #ifdef ROBUST_INTEGRAL_CHECK
    if (v002>zero2 &&
        abs(r0000+u00)>zero &&
        abs(r0100+u01)>zero &&
        abs(r1000+u10)>zero &&
        abs(r1100+u11)>zero    ){
    #else
    if (v002>zero2){
    #endif
        temp += +v00*( +u00v00*(-1 + 2*log(r0000+u00))
                       -u01v00*(-1 + 2*log(r0100+u01))
                       -u10v00*(-1 + 2*log(r1000+u10))
                       +u11v00*(-1 + 2*log(r1100+u11)) );
    }

    #ifdef ROBUST_INTEGRAL_CHECK
    if (v012>zero2 &&
        abs(r0001+u00)>zero &&
        abs(r0101+u01)>zero &&
        abs(r1001+u10)>zero &&
        abs(r1101+u11)>zero    ){
    #else
    if (v012>zero2){
    #endif
        temp += -v01*( +u00v01*(-1 + 2*log(r0001+u00))
                       -u01v01*(-1 + 2*log(r0101+u01))
                       -u10v01*(-1 + 2*log(r1001+u10))
                       +u11v01*(-1 + 2*log(r1101+u11)) );
    }

    #ifdef ROBUST_INTEGRAL_CHECK
    if (v102>zero2 &&
        abs(r0010+u00)>zero &&
        abs(r0110+u01)>zero &&
        abs(r1010+u10)>zero &&
        abs(r1110+u11)>zero    ){
    #else
    if (v102>zero2){
    #endif
        temp += -v10*( +u00v10*(-1 + 2*log(r0010+u00))
                       -u01v10*(-1 + 2*log(r0110+u01))
                       -u10v10*(-1 + 2*log(r1010+u10))
                       +u11v10*(-1 + 2*log(r1110+u11)) );
    }

    #ifdef ROBUST_INTEGRAL_CHECK
    if (v112>zero2 &&
        abs(r0011+u00)>zero &&
        abs(r0111+u01)>zero &&
        abs(r1011+u10)>zero &&
        abs(r1111+u11)>zero    ){
    #else
    if (v112>zero2){
    #endif
        temp += +v11*( +u00v11*(-1 + 2*log(r0011+u00))
                       -u01v11*(-1 + 2*log(r0111+u01))
                       -u10v11*(-1 + 2*log(r1011+u10))
                       +u11v11*(-1 + 2*log(r1111+u11)) );
    }
    value += 3*temp;

    #ifdef ROBUST_INTEGRAL_CHECK
    if ( z2 > zero2 &&
         abs(r0000-v00)>zero &&
         abs(r0100-v00)>zero &&
         abs(r1000-v00)>zero &&
         abs(r1100-v00)>zero &&

         abs(r0001-v01)>zero &&
         abs(r0101-v01)>zero &&
         abs(r1001-v01)>zero &&
         abs(r1101-v01)>zero &&

         abs(r0010-v10)>zero &&
         abs(r0110-v10)>zero &&
         abs(r1010-v10)>zero &&
         abs(r1110-v10)>zero &&

         abs(r0011-v11)>zero &&
         abs(r0111-v11)>zero &&
         abs(r1011-v11)>zero &&
         abs(r1111-v11)>zero &&

         abs(r0000-u00)>zero &&
         abs(r0001-u00)>zero &&
         abs(r0010-u00)>zero &&
         abs(r0011-u00)>zero &&

         abs(r0100-u01)>zero &&
         abs(r0101-u01)>zero &&
         abs(r0110-u01)>zero &&
         abs(r0111-u01)>zero &&

         abs(r1000-u10)>zero &&
         abs(r1001-u10)>zero &&
         abs(r1010-u10)>zero &&
         abs(r1011-u10)>zero &&

         abs(r1100-u11)>zero &&
         abs(r1101-u11)>zero &&
         abs(r1110-u11)>zero &&
         abs(r1111-u11)>zero     ){
    #else
    if ( z2 > zero2){
    #endif


        value +=
                4*(
                    -3*z*(
                        +u00v00*atan(u00v00/(r0000*z))
                        -u00v01*atan(u00v01/(r0001*z))
                        -u00v10*atan(u00v10/(r0010*z))
                        +u00v11*atan(u00v11/(r0011*z))

                        -u01v00*atan(u01v00/(r0100*z))
                        +u01v01*atan(u01v01/(r0101*z))
                        +u01v10*atan(u01v10/(r0110*z))
                        -u01v11*atan(u01v11/(r0111*z))

                        -u10v00*atan(u10v00/(r1000*z))
                        +u10v01*atan(u10v01/(r1001*z))
                        +u10v10*atan(u10v10/(r1010*z))
                        -u10v11*atan(u10v11/(r1011*z))

                        +u11v00*atan(u11v00/(r1100*z))
                        -u11v01*atan(u11v01/(r1101*z))
                        -u11v10*atan(u11v10/(r1110*z))
                        +u11v11*atan(u11v11/(r1111*z)) )


                    +z2*(
                        1.5* (
                            +v00*(
                                +log(r0000-v00)
                                -log(r0100-v00)
                                -log(r1000-v00)
                                +log(r1100-v00)	)
                            -v01*(
                                +log(r0001-v01)
                                -log(r0101-v01)
                                -log(r1001-v01)
                                +log(r1101-v01)	)
                            -v10*(
                                +log(r0010-v10)
                                -log(r0110-v10)
                                -log(r1010-v10)
                                +log(r1110-v10)	)
                            +v11*(
                                +log(r0011-v11)
                                -log(r0111-v11)
                                -log(r1011-v11)
                                +log(r1111-v11)	)

                            +u00*(
                                +log(r0000-u00)
                                -log(r0001-u00)
                                -log(r0010-u00)
                                +log(r0011-u00) )

                            -u01*(
                                +log(r0100-u01)
                                -log(r0101-u01)
                                -log(r0110-u01)
                                +log(r0111-u01) )

                            -u10*(
                                +log(r1000-u10)
                                -log(r1001-u10)
                                -log(r1010-u10)
                                +log(r1011-u10) )

                            +u11*(
                                +log(r1100-u11)
                                -log(r1101-u11)
                                -log(r1110-u11)
                                +log(r1111-u11) ) )

                        +r0000
                        -r0001
                        -r0010
                        +r0011

                        -r0100
                        +r0101
                        +r0110
                        -r0111

                        -r1000
                        +r1001
                        +r1010
                        -r1011

                        +r1100
                        -r1101
                        -r1110
                        +r1111
                        ));

    }
	return value/12;
}


float int_xyyz(float a, float b, float ly, float lz, float x, float y, float z){
    float xp[nBit];
    float yp[nBit];
    float zp[nBit];

    float xx[nBit];
    float yy[nBit];
    float zz[nBit];

    xp[MIN] = -a/2;
    xp[MAX] =  a/2;
    yp[MIN] = -b/2;
    yp[MAX] =  b/2;
    zp[MIN] =  0;
    zp[MAX] =  0;

	xp[CENTER] = 0;
	yp[CENTER] = 0;
	zp[CENTER] = 0;

	xp[LENGTH] = a;
	yp[LENGTH] = b;
	zp[LENGTH] = 0;


    xx[MIN] =  x;
    xx[MAX] =  x;
    yy[MIN] = -ly/2 +y;
    yy[MAX] =  ly/2 +y;
    zz[MIN] = -lz/2 +z;
    zz[MAX] =  lz/2 +z;

	xx[CENTER] = x;
	yy[CENTER] = y;
	zz[CENTER] = z;

	xx[LENGTH] = 0;
	yy[LENGTH] = ly;
	zz[LENGTH] = lz;

    float* coord1[nDim][nBit];
    float* coord2[nDim][nBit];

    *coord1[X] = xp;
    *coord1[Y] = yp;
    *coord1[Z] = zp;

    *coord2[X] = xx;
    *coord2[Y] = yy;
    *coord2[Z] = zz;

    return int_xyyz(coord1, coord2);
}


float int_xyyz(float* p1[3][4], float* p2[3][4]){

    #ifdef CAPLET_ATAN_LOG_INT_XYYZ
    using caplet::atan;
    using caplet::log;
    #else
    using std::atan;
    using std::log;
    #endif

	float* xp;
	float* yp;
	float  zp;

	float  xx;
	float* yy;

	float* zz;

    //* Make xp.length > zz.length and yp.length > yy.length
	if ( (*p1[X])[LENGTH] > (*p2[Z])[LENGTH] ){
		xp = *p1[X];
		zp = (*p1[Z])[CENTER];
		xx = (*p2[X])[CENTER];
		zz = *p2[Z];
	}else{
		xp = *p2[Z];
		zp = (*p2[X])[CENTER];
		xx = (*p1[Z])[CENTER];
		zz = *p1[X];
	}
	if ( (*p1[Y])[LENGTH] > (*p2[Y])[LENGTH] ){
		yp = *p1[Y];
		yy = *p2[Y];
	}else{
		yp = *p2[Y];
		yy = *p1[Y];
	}

	float a= xp[LENGTH];
	float b= yp[LENGTH];
    float ly=yy[LENGTH];
	float lz=zz[LENGTH];
	float x = abs( xp[CENTER] - xx );
	float y = abs( yp[CENTER] - yy[CENTER] );
	float z = abs( zp - zz[CENTER] );


    //* Approximation when panels are far seperated
	if ( b>a ){ // type1 criterion
		float ba = b/a;
		float lyb = ly/b;

		if ( ba > 15 ){
			if (
					x > ( ( lyb < 0.3 )? ( 0.4*b ) : (-0.25333+2.6695*lyb )*b )
						||
					y > ( 0.76324 + ( 1.7926 + 1.0448*lyb )*lyb )*b
						||
					z > ( 0.193145 + ( -0.41222 +2.5531*lyb )*lyb )*b
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}else if ( ba > 8 ){
			if (
					x > ( ( lyb < 0.3 )? ( 0.4*b ) : (-0.25333+2.6695*lyb )*b )
						||
					y > ( 0.76324 + ( 1.7926 + 1.0448*lyb )*lyb )*b
						||
					z > ( (0.00041439-1.8114e-05*ba+0.0005263*lyb)/(0.0014429-1.3985e-05*ba-0.00097839*lyb) +0.1 )*b
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}else if ( ba > 4 ){
			if (
					x > ( ( lyb < 0.3 )? ( 0.6*b ) : ( 0.1101 + 2.4728*lyb )*b)
						||
					y > ( 1.01025+ ( 0.83062 +1.8885*lyb)*lyb )*b
						||
					z > ( 2.1753 + (-0.38629 + 0.022515 *ba )*ba + (-4.0881 + 3.6737 *lyb +0.34002*ba )*lyb  )*b
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}else{ // ba > 1
			if (
					x > ( 3.7136+ (-1.6991 + 0.23106*ba )*ba + (-0.3211 + 1.0531*lyb +0.29667*ba )*lyb )*b
						||
					y > ( 3.4775 + ( -1.3063 + 0.15819*ba )*ba + ( -4.6679 + 3.3051 *lyb +1.2226*ba )*lyb )*b
						||
					z > ( 2.5 )*b
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}
	}else{ // type2 criterion
		float ab = a/b;
		float lza = lz/a;

		if ( ab > 8 ){
			if (
					x > ( 0.72357+1.7797*lza  )*a
						||
					y > ( 0.35665 + (-0.0061875 + 8.5796e-05*ab )*ab + ( 1.7517 +0.37208*lza + 0.0055581*ab )*lza )*a
						||
					z > ( ( 0.00030669+0.0068952*lza)/( 0.003384-0.0011036*lza ) + 0.1 )*a
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}else if ( ab > 4 ){
			if (
					x > ( 0.84223 - 0.031616*ab + 1.6364*lza )*a
						||
					y > ( 1.5194 +(-0.27035+0.015886*ab)*ab + (-1.0227+1.8757*lza+0.17602*ab)*lza )*a
						||
					z > ( 1.4212 + (-0.31988 + 0.022624*ab )*ab + ( 0.1677 + 2.1045 *lza + 0.13869 *ab )*lza  )*a
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}else if ( ab > 1.5 ){
			if (
					x > ( 2.6614 + ( -0.90581 + 0.11298*ab )*ab + ( 0.14021 + 0.93265*lza + 0.16342*ab )*lza )*a
						||
					y > ( 2.33 )*a
						||
					z > ( 2.574 + ( -0.80825 + 0.06658 *ab )*ab + (-5.059 +3.6317*lza +1.277*ab )*lza )*a
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}else{ // ab > 1
			if (
					x > ( 2.8 )*a
						||
					y > ( 3 )*a
						||
					z > ( 2.5 )*a
			){
				return int_xy( a, b, x, y, z, a*b )*lz*ly;
			}
		}
	}


    //* Analytical integral
    //	double value = 0.0;
    //	for (int m=0; m<2; m++){
    //	    double u = xp[m]-xx;
    //	    double u_2 = u*u;
    //	    for (int p=0; p<2; p++){
    //	        for (int q=0; q<2; q++){
    //	            double v = yp[p]-yy[q];
    //	            double v_2 = v*v;
    //	            for (int h=0; h<2; h++){
    //	                double z = zz[h]-zp;
    //	                double z_2 = z*z;
    //	                double uz = u*z;
    //
    //	                double r = sqrt(u_2+v_2+z_2);
    //	                double h11 =  2*uz*r ;
    //
    //	                if (v_2 > zero2){
    //	                     h11 += v_2*( +v*atan(u*z/v/r) );
    //	                }
    //
    //	                if (u_2 > zero2){
    //	                     h11 +=
    //	                         +u*( u_2 -3*v_2 ) * log(z + r)
    //	                         +3*v* ( u_2*atan(v*z/u/r)
    //	                                +2*uz*(1-log(v+r)) );
    //	                }
    //
    //	                if (z_2 > zero2){
    //	                     h11 +=
    //	                         +z*( z_2 -3*v_2 ) * log(u + r)
    //	                         +3*v*z_2*atan(u*v/z/r);
    //	                }
    //
    //	                value += (1-2*((m+p+q+h)%2)) *h11;
    //	    		}
    //	    	}
    //	    }
    //	}

    //* Expanded analytical integral
	double u1 = xp[0]-xx; double u1_2 = u1*u1;
	double u2 = xp[1]-xx; double u2_2 = u2*u2;

	double v11 = yp[0]-yy[0]; double v11_2 = v11*v11;
	double v12 = yp[0]-yy[1]; double v12_2 = v12*v12;
	double v21 = yp[1]-yy[0]; double v21_2 = v21*v21;
	double v22 = yp[1]-yy[1]; double v22_2 = v22*v22;

	double z1 = zz[0]-zp; double z1_2 = z1*z1;
	double z2 = zz[1]-zp; double z2_2 = z2*z2;

	double u1_z1 = u1*z1;
	double u1_z2 = u1*z2;
	double u2_z1 = u2*z1;
	double u2_z2 = u2*z2;

	double r1111 = sqrt( u1_2 + v11_2 + z1_2 );
	double r1112 = sqrt( u1_2 + v11_2 + z2_2 );
	double r1121 = sqrt( u1_2 + v12_2 + z1_2 );
	double r1122 = sqrt( u1_2 + v12_2 + z2_2 );

	double r1211 = sqrt( u1_2 + v21_2 + z1_2 );
	double r1212 = sqrt( u1_2 + v21_2 + z2_2 );
	double r1221 = sqrt( u1_2 + v22_2 + z1_2 );
	double r1222 = sqrt( u1_2 + v22_2 + z2_2 );

	double r2111 = sqrt( u2_2 + v11_2 + z1_2 );
	double r2112 = sqrt( u2_2 + v11_2 + z2_2 );
	double r2121 = sqrt( u2_2 + v12_2 + z1_2 );
	double r2122 = sqrt( u2_2 + v12_2 + z2_2 );

	double r2211 = sqrt( u2_2 + v21_2 + z1_2 );
	double r2212 = sqrt( u2_2 + v21_2 + z2_2 );
	double r2221 = sqrt( u2_2 + v22_2 + z1_2 );
	double r2222 = sqrt( u2_2 + v22_2 + z2_2 );


	double value =
			+2*(
					+u1_z1* ( r1111 - r1121 - r1211 + r1221 )
					-u1_z2* ( r1112 - r1122 - r1212 + r1222 )
					-u2_z1* ( r2111 - r2121 - r2211 + r2221 )
					+u2_z2* ( r2112 - r2122 - r2212 + r2222 )
				);

	if ( v11_2 > zero2 ){
		value +=
				+v11_2*v11*(
						+atan( u1*z1/(v11*r1111))
						-atan( u1*z2/(v11*r1112))
						-atan( u2*z1/(v11*r2111))
						+atan( u2*z2/(v11*r2112))
				);
	}
	if ( v12_2 > zero2 ){
		value -=
				+v12_2*v12*(
						+atan( u1*z1/(v12*r1121))
						-atan( u1*z2/(v12*r1122))
						-atan( u2*z1/(v12*r2121))
						+atan( u2*z2/(v12*r2122))
				);
	}
	if ( v21_2 > zero2 ){
		value -=
				+v21_2*v21*(
						+atan( u1*z1/(v21*r1211))
						-atan( u1*z2/(v21*r1212))
						-atan( u2*z1/(v21*r2211))
						+atan( u2*z2/(v21*r2212))
				);
	}
	if ( v22_2 > zero2 ){
		value +=
				+v22_2*v22*(
						+atan( u1*z1/(v22*r1221))
						-atan( u1*z2/(v22*r1222))
						-atan( u2*z1/(v22*r2221))
						+atan( u2*z2/(v22*r2222))
				);
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if (u1_2 > zero2 &&
        abs(z1+r1111)>zero  && abs(z2+r1112)>zero  &&
        abs(z1+r1121)>zero  && abs(z2+r1122)>zero  &&
        abs(z1+r1211)>zero  && abs(z2+r1212)>zero  &&
        abs(z1+r1221)>zero  && abs(z2+r1222)>zero  &&
		abs(v11+r1111)>zero && abs(v11+r1112)>zero && 
		abs(v12+r1121)>zero && abs(v12+r1122)>zero && 
		abs(v21+r1211)>zero && abs(v21+r1212)>zero && 
		abs(v22+r1221)>zero && abs(v22+r1222)>zero    ){
	#else
	if (u1_2 > zero2){
	#endif
		value += u1*(
				+( u1_2 - 3*v11_2 )*log( (z1+r1111)/(z2+r1112) )
				-( u1_2 - 3*v12_2 )*log( (z1+r1121)/(z2+r1122) )
				-( u1_2 - 3*v21_2 )*log( (z1+r1211)/(z2+r1212) )
				+( u1_2 - 3*v22_2 )*log( (z1+r1221)/(z2+r1222) )
				)
				+3*(
					+v11*(
						u1_2* ( atan( v11*z1/(u1*r1111) ) - atan( v11*z2/(u1*r1112) ) )
					    +2* ( u1_z1* ( 1-log(v11+r1111) ) - u1_z2* ( 1-log(v11+r1112) ) )
					)
					-v12*(
						u1_2* ( atan( v12*z1/(u1*r1121) ) - atan( v12*z2/(u1*r1122) ) )
					    +2* ( u1_z1* ( 1-log(v12+r1121) ) - u1_z2* ( 1-log(v12+r1122) ) )
					)
					-v21*(
						u1_2* ( atan( v21*z1/(u1*r1211) ) - atan( v21*z2/(u1*r1212) ) )
					    +2* ( u1_z1* ( 1-log(v21+r1211) ) - u1_z2* ( 1-log(v21+r1212) ) )
					)
					+v22*(
						u1_2* ( atan( v22*z1/(u1*r1221) ) - atan( v22*z2/(u1*r1222) ) )
					    +2* ( u1_z1* ( 1-log(v22+r1221) ) - u1_z2* ( 1-log(v22+r1222) ) )
					)
				);
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if (u2_2 > zero2 &&
        abs(z1+r2111)>zero  && abs(z2+r2112)>zero  &&
        abs(z1+r2121)>zero  && abs(z2+r2122)>zero  &&
        abs(z1+r2211)>zero  && abs(z2+r2212)>zero  &&
        abs(z1+r2221)>zero  && abs(z2+r2222)>zero  &&
		abs(v11+r2111)>zero && abs(v11+r2112)>zero && 
		abs(v12+r2121)>zero && abs(v12+r2122)>zero && 
		abs(v21+r2211)>zero && abs(v21+r2212)>zero &&
		abs(v22+r2221)>zero && abs(v22+r2222)>zero    ){
	#else
	if (u2_2 > zero2){
	#endif
		value -= u2*(
				+( u2_2 - 3*v11_2 )*log( (z1+r2111)/(z2+r2112) )
				-( u2_2 - 3*v12_2 )*log( (z1+r2121)/(z2+r2122) )
				-( u2_2 - 3*v21_2 )*log( (z1+r2211)/(z2+r2212) )
				+( u2_2 - 3*v22_2 )*log( (z1+r2221)/(z2+r2222) )
				)
				+3*(
					+v11*(
						u2_2* ( atan( v11*z1/(u2*r2111) ) - atan( v11*z2/(u2*r2112) ) )
					    +2* ( u2_z1* ( 1-log(v11+r2111) ) - u2_z2* ( 1-log(v11+r2112) ) )
					)
					-v12*(
						u2_2* ( atan( v12*z1/(u2*r2121) ) - atan( v12*z2/(u2*r2122) ) )
					    +2* ( u2_z1* ( 1-log(v12+r2121) ) - u2_z2* ( 1-log(v12+r2122) ) )
					)
					-v21*(
						u2_2* ( atan( v21*z1/(u2*r2211) ) - atan( v21*z2/(u2*r2212) ) )
					    +2* ( u2_z1* ( 1-log(v21+r2211) ) - u2_z2* ( 1-log(v21+r2212) ) )
					)
					+v22*(
						u2_2* ( atan( v22*z1/(u2*r2221) ) - atan( v22*z2/(u2*r2222) ) )
					    +2* ( u2_z1* ( 1-log(v22+r2221) ) - u2_z2* ( 1-log(v22+r2222) ) )
					)
				);
	}


	#ifdef ROBUST_INTEGRAL_CHECK
	if (z1_2 > zero2 &&
		abs(u1 + r1111)>zero && abs(u2 + r2111)>zero &&
		abs(u1 + r1121)>zero && abs(u2 + r2121)>zero && 
		abs(u1 + r1211)>zero && abs(u2 + r2211)>zero &&
		abs(u1 + r1221)>zero && abs(u2 + r2221)>zero    ){
	#else
	if (z1_2 > zero2){
	#endif
		value += z1*(
				+( z1_2 - 3*v11_2 )*log( (u1 + r1111)/(u2 + r2111) )
				-( z1_2 - 3*v12_2 )*log( (u1 + r1121)/(u2 + r2121) )
				-( z1_2 - 3*v21_2 )*log( (u1 + r1211)/(u2 + r2211) )
				+( z1_2 - 3*v22_2 )*log( (u1 + r1221)/(u2 + r2221) )
				)
				+3*(
					+v11*z1_2* ( atan( v11*u1/(z1*r1111) ) - atan( v11*u2/(z1*r2111) ) )
					-v12*z1_2* ( atan( v12*u1/(z1*r1121) ) - atan( v12*u2/(z1*r2121) ) )
					-v21*z1_2* ( atan( v21*u1/(z1*r1211) ) - atan( v21*u2/(z1*r2211) ) )
					+v22*z1_2* ( atan( v22*u1/(z1*r1221) ) - atan( v22*u2/(z1*r2221) ) )
				);
	}


	#ifdef ROBUST_INTEGRAL_CHECK
	if (z2_2 > zero2 &&
		abs(u1 + r1112)>zero && abs(u2 + r2112)>zero &&
		abs(u1 + r1122)>zero && abs(u2 + r2122)>zero && 
		abs(u1 + r1212)>zero && abs(u2 + r2212)>zero && 
		abs(u1 + r1222)>zero && abs(u2 + r2222)>zero 	){
	#else
	if (z2_2 > zero2){
	#endif
		value -= z2*(
				+( z2_2 - 3*v11_2 )*log( (u1 + r1112)/(u2 + r2112) )
				-( z2_2 - 3*v12_2 )*log( (u1 + r1122)/(u2 + r2122) )
				-( z2_2 - 3*v21_2 )*log( (u1 + r1212)/(u2 + r2212) )
				+( z2_2 - 3*v22_2 )*log( (u1 + r1222)/(u2 + r2222) )
				)
				+3*(
					+v11*z2_2* ( atan( v11*u1/(z2*r1112) ) - atan( v11*u2/(z2*r2112) ) )
					-v12*z2_2* ( atan( v12*u1/(z2*r1122) ) - atan( v12*u2/(z2*r2122) ) )
					-v21*z2_2* ( atan( v21*u1/(z2*r1212) ) - atan( v21*u2/(z2*r2212) ) )
					+v22*z2_2* ( atan( v22*u1/(z2*r1222) ) - atan( v22*u2/(z2*r2222) ) )
				);
	}

	return value /6;
}


float  int_xyzx(float* p1[nDim][nBit], float* p2[nDim][nBit]){
    float * p1Mirrored[nDim][nBit];
    float * p2Mirrored[nDim][nBit];
    *p1Mirrored[X] = *p1[Y];
    *p1Mirrored[Y] = *p1[X];
    *p1Mirrored[Z] = *p1[Z];
    *p2Mirrored[X] = *p2[Y];
    *p2Mirrored[Y] = *p2[X];
    *p2Mirrored[Z] = *p2[Z];

    return int_xyyz(p1Mirrored, p2Mirrored);
}


float int_xyy(float a, float b, float ly, float x, float y, float z){
    //**** ASSERT
    assert( a>0 );
	assert( b>0 );
	assert( ly>0 );
	assert( x>=0 );
	assert( y>=0 );
	assert( z>=0 );
    //****

    #ifdef CAPLET_ATAN_LOG_INT_XYY
    using caplet::atan;
    using caplet::log;
    #else
    using std::atan;
    using std::log;
    #endif

    //* Approximation when seperation is far
    if ( b >= a ){
        float lyb 	= ly/b;

		if ( lyb > 0.3 ){
			if (
					x > ( -0.03537+2.3839*lyb )*b
						||
					y > ( 0.19506 +3.1565*lyb )*b
						||
					z > ( -0.22164+2.4932*lyb )*b
			){
				return int_xy(a,b,x,y,z,a*b)*ly;
			}
		}else{ // lyb > 0
			return int_xy(a,b,x,y,z,a*b)*ly;
		}
	}else{ // b < a
		float ab 	= a/b;
		float lyb	= ly/b;

		if ( lyb > 0.3 ){
			if ( ab > 8 ){
				if (
						x > ( 0.6 )*a
							||
						y > ( (-8.6945e-05 + 4.6425e-06*ab + 0.00048786*lyb)/(-0.0016195 + 0.00038443*ab + 0.0001899 *lyb) )*a
						 	||
						z > ( (-1.8871e-05 - 2.7441e-07*ab + 9.6839e-05*lyb)/(-1.8457e-05+ 4.5722e-05*ab + 8.3233e-06*lyb) )*a
				){
					return int_xy(b,a,y,x,z,a*b)*ly;
				}
			}else{ // ab > 1
				if (
				        x > ( (-0.00015075 + 0.00012366*ab + 0.00046419*lyb)/(-8.9475e-05 + 0.00026222*ab + 1.4625e-05*lyb) )*a
				        	||
				        y > ( ( 3.0636e-05 - 1.2516e-05*ab + 0.00043222*lyb)/(-3.4346e-05 + 0.0001804 *ab - 1.944e-05 *lyb) )*a
				        	||
				        z > ( (-1.8871e-05 - 2.7441e-07*ab + 9.6839e-05*lyb)/(-1.8457e-05 + 4.5722e-05*ab + 8.3233e-06*lyb) )*a
				){
					return int_xy(b,a,y,x,z,a*b)*ly;
				}
			}
		}else if ( lyb > 0.2 ){
			if ( ab > 8 ){
				if ( y > ( (-8.6945e-05 + 4.6425e-06*ab + 0.00048786*lyb)/(-0.0016195 + 0.00038443*ab + 0.0001899*lyb) )*a ){
					return int_xy(b,a,y,x,z,a*b)*ly;
				}
			}else{ // ab > 1
				if ( y > ( (3.0636e-05-1.2516e-05*ab + 0.00043222*lyb)/(-3.4346e-05 + 0.0001804*ab -1.944e-05*lyb) )*a ){
					return int_xy(b,a,y,x,z,a*b)*ly;
				}
			}
		}else{ // lyb > 0
			return int_xy(b,a,y,x,z,a*b)*ly;
		}
	}

    //* Analytical integral
	float xp[] = {-a/2,  a/2};
	float yp[] = {-b/2,  b/2};
	float xx   = x;
	float yy[] = {-ly/2+y, ly/2+y};
	float z_2 = z*z;

	float u1 = xp[0]-xx;		float u1_2 = u1*u1;
	float u2 = xp[1]-xx;		float u2_2 = u2*u2;

	float v11 = yp[0]-yy[0];	float v11_2 = v11*v11;
	float v12 = yp[0]-yy[1];	float v12_2 = v12*v12;
	float v21 = yp[1]-yy[0];	float v21_2 = v21*v21;
	float v22 = yp[1]-yy[1];	float v22_2 = v22*v22;

	float r111 = sqrt( u1_2 + v11_2 + z_2 );
	float r112 = sqrt( u1_2 + v12_2 + z_2 );
	float r121 = sqrt( u1_2 + v21_2 + z_2 );
	float r122 = sqrt( u1_2 + v22_2 + z_2 );

	float r211 = sqrt( u2_2 + v11_2 + z_2 );
	float r212 = sqrt( u2_2 + v12_2 + z_2 );
	float r221 = sqrt( u2_2 + v21_2 + z_2 );
	float r222 = sqrt( u2_2 + v22_2 + z_2 );


	float val = 0;



    val = -( +u1*(r111-r112-r121+r122)
             -u2*(r211-r212-r221+r222) );

	if ( z_2 > zero2 ){
		val += -2*z*(
					+v11*(  atan(u1*v11/(z*r111)) - atan(u2*v11/(z*r211)) )
					-v12*(  atan(u1*v12/(z*r112)) - atan(u2*v12/(z*r212)) )
					-v21*(  atan(u1*v21/(z*r121)) - atan(u2*v21/(z*r221)) )
					+v22*(  atan(u1*v22/(z*r122)) - atan(u2*v22/(z*r222)) )
				);
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if ( abs(v11_2-z_2) > zero2 && 
		 abs(u1+r111)>zero && 
		 abs(u2+r211)>zero    ){
	#else
	if ( abs(v11_2-z_2) > zero2 ){
	#endif
		val += ( v11_2 - z_2 )*( log( (u1+r111)/(u2+r211) ) );
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if ( abs(v12_2-z_2) > zero2 &&
		 abs(u1+r112)>zero &&
		 abs(u2+r212)>zero    ){
	#else
	if ( abs(v12_2-z_2) > zero2 ){
	#endif
		val -= ( v12_2 - z_2 )*( log( (u1+r112)/(u2+r212) ) );
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if ( abs(v21_2-z_2) > zero2 &&
		 abs(u1+r121)>zero &&
		 abs(u2+r221)>zero    ){
	#else
	if ( abs(v21_2-z_2) > zero2 ){
	#endif
		val -= ( v21_2 - z_2 )*( log( (u1+r121)/(u2+r221) ) );
	}

	#ifdef ROBUST_INTEGRAL_CHECK	
	if ( abs(v22_2-z_2) > zero2 &&
		 abs(u1+r122)>zero &&
		 abs(u2+r222)>zero    ){
	#else
	if ( abs(v22_2-z_2) > zero2 ){
	#endif
		val += ( v22_2 - z_2 )*( log( (u1+r122)/(u2+r222) ) );
	}

    #ifdef ROBUST_INTEGRAL_CHECK
    if ( u1_2 > zero2 &&
         abs(v11+r111)>zero &&
         abs(v12+r112)>zero &&
         abs(v21+r121)>zero &&
         abs(v22+r122)>zero    ){
    #else
    if ( u1_2 > zero2 ){
    #endif
        val += 2*u1* (
					+ v11*log(v11+r111)
					- v12*log(v12+r112)
					- v21*log(v21+r121)
					+ v22*log(v22+r122)
				);
	}

    #ifdef ROBUST_INTEGRAL_CHECK
    if ( u2_2 > zero2 &&
         abs(v11+r211)>zero &&
         abs(v12+r212)>zero &&
         abs(v21+r221)>zero &&
         abs(v22+r222)>zero    ){
    #else
    if ( u2_2 > zero2 ){
    #endif
        val -= 2*u2* (
					+ v11*log(v11+r211)
					- v12*log(v12+r212)
					- v21*log(v21+r221)
					+ v22*log(v22+r222)
				);
	}

	return val/2;
}


float int_xyz(float a, float b, float lz, float x, float y, float z){
    //**** ASSERT
	assert( a>0 );
	assert( b>0 );
	assert( lz>0 );
	assert( x>=0 );
	assert( y>=0 );
	assert( z>=0 );
    //****

    #ifdef CAPLET_ATAN_LOG_INT_XYZ
    using caplet::atan;
    using caplet::log;
    #else
    using std::atan;
    using std::log;
    #endif


    //* Rotate to make sure b > a > lz
    //  in order to make effective approximation when possible
	float temp;
	if( b > a ){
	    if( lz > a ){
	        if( b > lz ){ // b  > lz > a
	            temp = lz; lz = a; a = temp;
	            temp = z ; z  = x; x = temp;
	        }else{        // lz > b  > a
	            temp = a; a = b; b = lz; lz = temp;
	            temp = x; x = y; y = z ; z  = temp;
	        }
	    }
	}else{ // a > b
	    if( a > lz ){
	        if( b > lz ){  // a  > b  > lz
	            temp = a; a = b; b = temp;
	            temp = x; x = y; y = temp;
	        }else{        // a  > lz > b
	            temp = a; a = lz; lz=b; b = temp;
	            temp = x; x = z;  z =y; y = temp;
	        }
	    }else{  // lz > a  > b
	        temp = lz; lz = b; b = temp;
	        temp = z ; z  = y; y = temp;
	    }
	}
    // until here, b > a > lz


    //* Approximation formula for b > a > lz
	float ba = b/a;
	float lza = lz/a;

	// x-dir
	if ( x > ( (7.4772e-06 + (7.997e-08 + 6.4163e-09*ba )*ba + (1.8516e-05 + 1.2866e-05*lza -5.7301e-07*ba)*lza)/(-2.0543e-06 + 1.8656e-05*ba -9.1847e-07*lza) )*b ){
		return int_xy(a,b,x,y,z, a*b)*lz;
	}
	// y-dir
	if ( y > ( ( -4.9867e-06+9.531e-06*ba+1.8941e-05*lza )/(-6.3555e-06+1.8856e-05*ba-2.7403e-06*lza) +0.1 )*b ){
		return int_xy(a,b,x,y,z, a*b)*lz;
	}
	// z-dir
	if ( ba > 12 ){
		if ( z > 0.1*b ){
			return int_xy(a,b,x,y,z, a*b)*lz;
		}
	}else if ( ba > 3.5 ){
		if ( z > 0.6*b ){
			return int_xy(a,b,x,y,z, a*b)*lz;
		}
	}else{ // ba >= 1
		if ( z > 1.2*b ){
			return int_xy(a,b,x,y,z, a*b)*lz;
		}
	}

    //* Analytical integral
    //	float xp[] = {-a/2-x, a/2-x};
    //	float yp[] = {-b/2-y, b/2-y};
    //	float zz[] = {-lz/2-z, lz/2-z};
    //
    //	for ( int h=0; h<2; h++ ){
    //	      for ( int m = 0; m<2; m++ ){
    //	          for ( int p = 0; p<2; p++){
    //	              float u = xp[m]; float u_2 = u*u;
    //	              float v = yp[p]; float v_2 = v*v;
    //	              float z = zz[h]; float z_2 = z*z;
    //	              float r = sqrt(u*u+v*v+z*z);
    //
    //	              float temp = 0;
    //	              if ( z_2 > zero2 ){
    //	                  temp = temp - 0.5*z_2*atan(u*v/z/r);
    //	                  temp = temp + z*(v*log(u+r) + u*log(v+r));
    //	              }
    //	              if ( v_2 > zero2 ){
    //	                  temp = temp - 0.5*v_2*atan(u*z/v/r);
    //	                  temp = temp + u*v*log(z+r);
    //	              }
    //	              if  ( u_2 > zero2 ){
    //	                  temp = temp - 0.5*u_2*atan(v*z/u/r);
    //	              }
    //	              val = val + (1-2*((m+p+h+1)%2))*temp;
    //	          }
    //	      }
    //	}

    //* Expanded analytical integral
    float val = 0 ;

    float u1 = -a/2  -x;	float u1_2 = u1*u1;
    float u2 =  a/2  -x;	float u2_2 = u2*u2;
    float v1 = -b/2  -y;	float v1_2 = v1*v1;
    float v2 =  b/2  -y;	float v2_2 = v2*v2;
	float z1 = -lz/2 -z;	float z1_2 = z1*z1;
	float z2 =  lz/2 -z;	float z2_2 = z2*z2;

	float r111 = sqrt( u1_2 + v1_2 + z1_2 );
	float r112 = sqrt( u1_2 + v1_2 + z2_2 );
	float r121 = sqrt( u1_2 + v2_2 + z1_2 );
	float r122 = sqrt( u1_2 + v2_2 + z2_2 );
	float r211 = sqrt( u2_2 + v1_2 + z1_2 );
	float r212 = sqrt( u2_2 + v1_2 + z2_2 );
	float r221 = sqrt( u2_2 + v2_2 + z1_2 );
	float r222 = sqrt( u2_2 + v2_2 + z2_2 );


	#ifdef ROBUST_INTEGRAL_CHECK
	if ( z1_2 > zero2 &&
		 abs(v1+r111)>zero && 
		 abs(v2+r121)>zero && 
		 abs(v1+r211)>zero && 
		 abs(v2+r221)>zero &&
	  	 abs(u1+r111)>zero && 
	  	 abs(u2+r211)>zero && 
	  	 abs(u1+r121)>zero && 
	  	 abs(u2+r221)>zero     ){ 
	#else
	if ( z1_2 > zero2 ){
	#endif 
	    val += + 0.5* z1_2 *(
					  + atan( u1*v1/z1/r111 )
					  - atan( u1*v2/z1/r121 )
					  - atan( u2*v1/z1/r211 )
					  + atan( u2*v2/z1/r221 )
				)
				+z1*( -u1*log((v1+r111)/(v2+r121)) + u2*log((v1+r211)/(v2+r221))
					  -v1*log((u1+r111)/(u2+r211)) + v2*log((u1+r121)/(u2+r221))
				);
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if ( z2_2 > zero2 && 
         abs(v1+r112)>zero && 
         abs(v2+r122)>zero && 
         abs(v1+r212)>zero && 
         abs(v2+r222)>zero &&
         abs(u1+r112)>zero && 
         abs(u2+r212)>zero && 
         abs(u1+r122)>zero && 
         abs(u2+r222)>zero     ){		
	#else
	if ( z2_2 > zero2 ){
	#endif
	    val -= 0.5* z2_2 *(
					  + atan( u1*v1/z2/r112 )
					  - atan( u1*v2/z2/r122 )
					  - atan( u2*v1/z2/r212 )
					  + atan( u2*v2/z2/r222 )
	          )
        	  + z2*(- u1*log((v1+r112)/(v2+r122)) + u2*log((v1+r212)/(v2+r222))
        			- v1*log((u1+r112)/(u2+r212)) + v2*log((u1+r122)/(u2+r222))
        	  );
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if ( v1_2 > zero2 &&
		 abs(z1+r111)>zero && 
		 abs(z2+r112)>zero && 
		 abs(z1+r211)>zero && 
		 abs(z2+r212)>zero    ){
	#else
	if ( v1_2 > zero2 ){
	#endif
	    val += 0.5*v1_2*(
							+ atan(u1*z1/v1/r111)
							- atan(u1*z2/v1/r112)
							- atan(u2*z1/v1/r211)
							+ atan(u2*z2/v1/r212)
					)
	    		+v1*( -u1*log((z1+r111)/(z2+r112))
	    			  +u2*log((z1+r211)/(z2+r212))
        );
	}

	#ifdef ROBUST_INTEGRAL_CHECK
	if ( v2_2 > zero2 &&
		 abs(z1+r121)>zero && 
		 abs(z2+r122)>zero && 
		 abs(z1+r221)>zero && 
		 abs(z2+r222)>zero     ){
	#else
	if ( v2_2 > zero2 ){
	#endif	
	    val -= 0.5*v2_2*(
							+ atan(u1*z1/v2/r121)
							- atan(u1*z2/v2/r122)
							- atan(u2*z1/v2/r221)
							+ atan(u2*z2/v2/r222)
					)
	    	    +v2*( -u1*log((z1+r121)/(z2+r122))
                      +u2*log((z1+r221)/(z2+r222))
        );
	}

	if ( u1_2 > zero2 ){
	    val += 0.5*u1_2*(
							+ atan(v1*z1/u1/r111)
							- atan(v1*z2/u1/r112)
							- atan(v2*z1/u1/r121)
							+ atan(v2*z2/u1/r122)
						);
	}
	if ( u2_2 > zero2 ){
	    val -= 0.5*u2_2*(
							+ atan(v1*z1/u2/r211)
							- atan(v1*z2/u2/r212)
							- atan(v2*z1/u2/r221)
							+ atan(v2*z2/u2/r222)
						);
	}
	return val;
}



//***********************************************************************
//*
//* DOUBLE VERSION
//*
//*

double int_xy_d(float a, float b, float x, float y, float z, float area);
double calColD(float* p1[nDim][nBit], float* p2[nDim][nBit]){
    //* Test point : panel1.center-panel2.center
    //  Integration: panel2
    float a = (*p2[X])[LENGTH];
    float b = (*p2[Y])[LENGTH];

    float x = (*p1[X])[CENTER] - (*p2[X])[CENTER];
    float y = (*p1[Y])[CENTER] - (*p2[Y])[CENTER];
    float z = (*p1[Z])[CENTER] - (*p2[Z])[CENTER];

    return int_xy_d(a,b,x,y,z,a*b);
}
double int_xy_d(float a, float b, float x, float y, float z, float area){
    //**** ASSERT
    assert( a>0 );
    assert( b>0 );
    assert( x>=0 );
    assert( y>=0 );
    assert( z>=0 );
    assert( area>0 );
    //****

    using std::atan;
    using std::log;

    //* Compute aspect ratio
    double ba = b/a;

    if ( ba < 1 ){ // swap x and y
        ba = 1/ba;
        float temp;
        temp = b; b = a; a = temp;
        temp = y; y = x; x = temp;
    }

    //* Approximation when the evaluation point is far from the origin
    //  add gr*b guard ring to improve accuracy
    if ( ba < 8 ){
        if ( x > (4.545f+approximationGuardRingForCalColD)*b
                 ||
             y > (( -0.036938f+0.044672f*ba )/( -0.0051053f+0.0068095f*ba )+approximationGuardRingForCalColD) *b
                 ||
             z > (( -0.011618f+0.029056f*ba )/( -0.0037776f+0.0064776f*ba )+approximationGuardRingForCalColD) *b
        ){
            return area/sqrt(x*x+y*y+z*z);
        }
    }else{
        if ( x > (4.545f+approximationGuardRingForCalColD)*b || y > (6.545f+approximationGuardRingForCalColD)*b || z > (4.636f+approximationGuardRingForCalColD)*b ){
            return area/sqrt(x*x+y*y+z*z);
        }
    }

    //* Analytical integral
    double x1 = x-a/2;
    double x2 = x+a/2;
    double y1 = y-b/2;
    double y2 = y+b/2;

    double x1_2 = x1*x1;
    double x2_2 = x2*x2;
    double y1_2 = y1*y1;
    double y2_2 = y2*y2;

    double z_2  = z*z;

    double r11 = sqrt(x1_2 + y1_2 + z_2);
    double r12 = sqrt(x1_2 + y2_2 + z_2);
    double r21 = sqrt(x2_2 + y1_2 + z_2);
    double r22 = sqrt(x2_2 + y2_2 + z_2);

    #ifdef ROBUST_INTEGRAL_CHECK
    double p = ( ( x2_2>zero2 && abs(y2+r22)>zero && abs(y1+r21)>zero )? x2*log( (y2+r22)/(y1+r21) ):0)
              +( ( x1_2>zero2 && abs(y1+r11)>zero && abs(y2+r12)>zero )? x1*log( (y1+r11)/(y2+r12) ):0)
              +( ( y2_2>zero2 && abs(x2+r22)>zero && abs(x1+r12)>zero )? y2*log( (x2+r22)/(x1+r12) ):0)
              +( ( y1_2>zero2 && abs(x1+r11)>zero && abs(x2+r21)>zero )? y1*log( (x1+r11)/(x2+r21) ):0)
              +( (  z_2>zero2 )?  z*( atan(x2*y1/r21/z) + atan(x1*y2/r12/z)
                                  -atan(x2*y2/r22/z) - atan(x1*y1/r11/z) ):0);
  	#else
    double p = ((x2_2>zero2)? x2*log( (y2+r22)/(y1+r21) ):0)
              +((x1_2>zero2)? x1*log( (y1+r11)/(y2+r12) ):0)
              +((y2_2>zero2)? y2*log( (x2+r22)/(x1+r12) ):0)
              +((y1_2>zero2)? y1*log( (x1+r11)/(x2+r21) ):0)
              +(( z_2>zero2)?  z*( atan(x2*y1/r21/z) + atan(x1*y2/r12/z)
                                  -atan(x2*y2/r22/z) - atan(x1*y1/r11/z) ):0);
    #endif
    return p;
}


} // end of namespace caplet
