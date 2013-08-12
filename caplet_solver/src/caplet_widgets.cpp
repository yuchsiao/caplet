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

#include "caplet_widgets.h"

#include <iostream>
#include <iomanip>

namespace caplet{

template <class T>
void print_vector(T* x, int nx, const std::string name){
	if (name.empty()==false) {
		std::cout << name << std::endl;
	}
	for (int i=0; i<nx; i++){
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;
}
template void print_vector<float> (float*  x, int nx, const std::string name);
template void print_vector<double>(double* x, int nx, const std::string name);
template void print_vector<int>   (int*    x, int nx, const std::string name);

template <class T>
void print_matrix(T* A, int m, int n, const std::string name, char uplo, std::ostream &out=std::cout){
	if (name.empty()==false) {
		out << name << std::endl;
	}
	// Fortran order

	if ( uplo == 'a' ){
		for (int i=0; i<m; i++){
			for (int j=0; j<n; j++){
                out << std::setw(14) << std::setprecision(6)  << A[i+m*j];
			}
            out << std::endl;
		}
	}else if ( uplo == 'u' ){
		for (int i=0; i<m; i++){
			for (int j=0; j<i; j++){
                out << std::setw(14) << std::setprecision(6) << 0;
			}
			for (int j=i; j<n; j++){
                out << std::setw(14)  << std::setprecision(6) << A[i+m*j];
			}
            out << std::endl;
		}
	}else if ( uplo == 'l' ){
		for (int i=0; i<m; i++){
			for (int j=0; j<=i; j++){
                out << std::setw(14) << std::setprecision(6) << A[i+m*j];
			}
			for (int j=i+1; j<n; j++){
                out << std::setw(14) << std::setprecision(6)  << 0.0f;
			}
            out << std::endl;
		}
	}
}
template void print_matrix<float> (float*  x, int m, int n, const std::string name, char uplo, std::ostream& out=std::cout);
template void print_matrix<double>(double* x, int m, int n, const std::string name, char uplo, std::ostream& out=std::cout);
template void print_matrix<int>   (int*    x, int m, int n, const std::string name, char uplo, std::ostream& out=std::cout);



}
