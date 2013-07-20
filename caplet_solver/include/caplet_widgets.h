/*
Created: Jul 27, 2010
Author : Yu-Chung Hsiao
Email  : project.caplet@gmail.com
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


#ifndef CAPLET_WIDGETS_H_
#define CAPLET_WIDGETS_H_

#include <string>

#include <iostream>
#include <fstream>

//        p1  = [0.];
//        w1  = [2.];
//    case 2
//        p1  = [-0.5773502691896257, ];
//        w1  = [1., 1.];
//    case 3
//        p1  = [-0.7745966692414834, ];
//        w1  = [0.5555555555555553, ];
//    case 4
//        p1  = [-0.861136311594053, -0.3399810435848562, ];
//        w1  = [0.3478548451374539, 0.6521451548625462, ];
//    case 5
//        p1  = [-0.906179845938664, -0.5384693101056829, ];
//        w1  = [0.2369268850561887, 0.4786286704993665, ];
//    case 6
//        p1  = [-0.932469514203152, -0.6612093864662646, -0.2386191860831968, ];
//        w1  = [0.1713244923791709, 0.3607615730481379, 0.4679139345726913, ];
//    case 7
//        p1  = [-0.949107912342759, -0.7415311855993937, -0.4058451513773972, ];
//        w1  = [0.129484966168868, 0.2797053914892783, 0.3818300505051186, ];
//    case 8
//        p1  = [-0.960289856497537, -0.7966664774136262, -0.5255324099163289, -0.1834346424956498, ];
//        w1  = [0.1012285362903738, 0.2223810344533786, 0.3137066458778874, 0.3626837833783619, ];
//    case 9
//        p1  = [-0.968160239507626, -0.836031107326637, -0.6133714327005903, -0.3242534234038088, ];
//        w1  = [0.0812743883615759, 0.1806481606948543, 0.2606106964029356, 0.3123470770400029, ];
//    case 10
//        p1  = [-0.973906528517172, -0.865063366688984, -0.6794095682990246, -0.433395394129247, -0.1488743389816312, ];
//        w1  = [0.06667134430868681, 0.149451349150573, 0.2190863625159832, 0.2692667193099968, 0.2955242247147529, ];
//    case 11
//        p1  = [-0.97822865814604, -0.88706259976812, -0.7301520055740422, -0.5190961292068116, -0.2695431559523449, ];
//        w1  = [0.05566856711621584, 0.1255803694648743, 0.1862902109277404, 0.2331937645919927, 0.2628045445102466, ];
//    case 12
//        p1  = [-0.981560634246732, -0.904117256370452, -0.7699026741943177, -0.5873179542866143, -0.3678314989981804, -0.1252334085114688, ];
//        w1  = [0.04717533638647547, 0.1069393259953637, 0.1600783285433586, 0.2031674267230672, 0.2334925365383534, 0.2491470458134027, ];


//const float gauss_p[10][5]
namespace caplet{

template <class T>
void print_vector(T* x, int nx, const std::string name="vector");
template <class T>
void print_matrix(T* A, int n, int m, const std::string name="matrix", char uplo = 'a', std::ostream& out=std::cout);



template <class T>
double qgaus(T &func, const double a, const double b){
	static const double x[]={
			0.1488743389816312,
			0.4333953941292472,
			0.6794095682990244,
			0.8650633666889845,
			0.9739065285171717
	};
	static const double w[]={
			0.2955242247147529,
			0.2692667193099963,
			0.2190863625159821,
			0.1494513491505806,
			0.0666713443086881
	};
	double xm=0.5*(b+a);
	double xr=0.5*(b-a);
	double s=0;
	for (int j=0;j<5;j++) {
		double dx=xr*x[j];
		s += w[j]*(func(xm+dx)+func(xm-dx));
	}
	return s *= xr;
}

template <class T>
float gaussquad1(T &func, const float a, const float b){
	return func((b+a)/2)*(b-a);
}
template <class T>
float gaussquad2(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( func( xm+xr*0.5773502691896257f ) + func( xm-xr*0.5773502691896257f ) )*xr;
}
template <class T>
float gaussquad3(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.888888888888889f*func(xm) + 0.5555555555555553f*( func(xm+xr*0.7745966692414834f) + func(xm-xr*0.7745966692414834f) ) )*xr;
}
template <class T>
float gaussquad4(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.6521451548625462f*(func(xm+xr*0.3399810435848562f)+func(xm-xr*0.3399810435848562f)) +
			 0.3478548451374539f*(func(xm+xr*0.861136311594053f) +func(xm-xr*0.861136311594053f )))*xr;
}
template <class T>
float gaussquad5(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.5688888888888889f* func(xm)
			+0.4786286704993665f* ( func( xm+xr*0.5384693101056829f )+ func( xm-xr*0.5384693101056829f ) )
			+0.2369268850561887f* ( func( xm+xr*0.906179845938664f )+ func( xm-xr*0.906179845938664f ) ))*xr;
}

template <class T>
float gaussquad6(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.4679139345726913f* ( func( xm+xr*0.2386191860831968f )+ func( xm-xr*0.2386191860831968f ) )
			+0.3607615730481379f* ( func( xm+xr*0.6612093864662646f )+ func( xm-xr*0.6612093864662646f ) )
			+0.1713244923791709f* ( func( xm+xr*0.932469514203152f )+ func( xm-xr*0.932469514203152f ) ))*xr;
}
template <class T>
float gaussquad7(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.4179591836734694f*   func( xm )
			+0.3818300505051188f* ( func( xm+xr*0.4058451513773971f )+ func( xm-xr*0.4058451513773971f ) )
			+0.279705391489276f * ( func( xm+xr*0.7415311855993945f )+ func( xm-xr*0.7415311855993945f ) )
			+0.1294849661688697f* ( func( xm+xr*0.949107912342759f  )+ func( xm-xr*0.949107912342759f  ) ))*xr;
}
template <class T>
float gaussquad8(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.3626837833783619f* ( func( xm+xr*0.1834346424956498f )+ func( xm-xr*0.1834346424956498f ) )
			+0.3137066458778874f* ( func( xm+xr*0.5255324099163289f )+ func( xm-xr*0.5255324099163289f ) )
			+0.2223810344533786f* ( func( xm+xr*0.7966664774136262f )+ func( xm-xr*0.7966664774136262f ) )
			+0.1012285362903738f* ( func( xm+xr*0.960289856497537f  )+ func( xm-xr*0.960289856497537f  ) ))*xr;
}
template <class T>
float gaussquad9(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.3302393550012597f*   func( xm )
			+0.3123470770400025f* ( func( xm+xr*0.3242534234038088f )+ func( xm-xr*0.3242534234038088f ) )
			+0.2606106964029353f* ( func( xm+xr*0.6133714327005908f )+ func( xm-xr*0.6133714327005908f ) )
			+0.1806481606948577f* ( func( xm+xr*0.836031107326635f  )+ func( xm-xr*0.836031107326635f  ) )
			+0.0812743883615721f* ( func( xm+xr*0.968160239507627f  )+ func( xm-xr*0.968160239507627f  ) ))*xr;
}
template <class T>
float gaussquad10(T &func, const float a, const float b){
	float xm = (b+a)/2;
	float xr = (b-a)/2;
	return ( 0.2955242247147529f* ( func( xm+xr*0.1488743389816312f )+ func( xm-xr*0.1488743389816312f ) )
			+0.2692667193099968f* ( func( xm+xr*0.433395394129247f  )+ func( xm-xr*0.433395394129247f  ) )
			+0.2190863625159832f* ( func( xm+xr*0.6794095682990246f )+ func( xm-xr*0.6794095682990246f ) )
			+0.149451349150573f * ( func( xm+xr*0.865063366688984f  )+ func( xm-xr*0.865063366688984f  ) )
			+0.06667134430868681f*( func( xm+xr*0.973906528517172f  )+ func( xm-xr*0.973906528517172f  ) ))*xr;
}
template <class T>
float gaussquad(T &func, const float a, const float b, const int n){
	// n-point gaussquad integration of func between a and b
	switch(n){
	case 1: return gaussquad1(func, a,b);
	case 2: return gaussquad2(func, a,b);
	case 3: return gaussquad3(func, a,b);
	case 4: return gaussquad4(func, a,b);
	case 5: return gaussquad5(func, a,b);
	default:
		switch(n){
		case 6: return gaussquad6(func, a,b);
		case 7: return gaussquad7(func, a,b);
		case 8: return gaussquad8(func, a,b);
		case 9: return gaussquad9(func, a,b);
		case 10: return gaussquad10(func, a,b);
		default:
			std::cerr << "ERROR: out of available degrees in the gauss quadrature" << std::endl;
		}
	}
	return 0.0f;
}



}

#endif /* CAPLET_WIDGETS_H_ */
