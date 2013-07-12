/*
Created: Jul 27, 2010
Author : Yu-Chung Hsiao
Email  : yuchsiao@mit.edu
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

#ifndef CAPLET_BLAS_H_
#define CAPLET_BLAS_H_

namespace caplet{

extern "C"{
	/* Blas part */
	// double precision dot product:
	double ddot_(
			const int* 		n, 		// dimension
			const double* 	dx, 	// []vector x
			const int* 		incx, 	// index increment of each access of x
			const double* 	dy, 	// []vector y
			const int* 		incy	// index increment of each access of y
	);
	// single precision dot product
	float sdot_(
			const int* 		n, 		// dimension
			const float* 	dx, 	// []vector x
			const int* 		incx, 	// index increment of each access of x
			const float* 	dy, 	// []vector y
			const int* 		incy	// index increment of each access of y
	);

	// absolute value sum: asum <- ||x||_1
	float sasum_(
			const int*		n,		// i	dimension
			const float*	sx,		// i[]
			const int*		incx	// i
	);

	// return the maximum absolute value of x
    double damax_(
    		const int*		n,		// i
    		const double*	dx,		// i[] n-by-1 vector
    		const int*		incx	// i   index increment of each access of x
    );

    // y := a * x + y
    void daxpy_(
    		const int*		n,		// i
    		const double*	da,		// i
    		const double*	dx,		// i[]
    		const int*		incx,	// i
    		double*			dy,		// io[]
    		const int*		incy	// i
    );
    void saxpy_(
    		const int*		n,		// i
    		const float*	da,		// i
    		const float*	dx,		// i[]
    		const int*		incx,	// i
    		float*			dy,		// io[]
    		const int*		incy	// i
    );

    // dy <- dx
    void dcopy_(
    		const int*		n,		// i
			const double*	dx,		// i[]
			const int*		incx,	// i
			double*			dy,		// o[]
			const int*		incy	// i
    );
    void scopy_(
    		const int*		n,		// i
    		const float*	sx,		// i[]
    		const int*		incx,	// i
    		float*			sy,		// o[]
    		const int*		incy	// i
    );

    // dx <- da dx
    void dscal_(
    		const int*		n,		// i
    		const double*	da,		// i
    		double*			dx,		// io[]
    		const int*		incx	// i
    );
    void sscal_(
    		const int*		n,		// i
    		const float*	sa,		// i
    		float*			sx,		// io[]
    		const int*		incx	// i
    );

    // A <- alpha x y' + A
    void dger_(
    		const int*		m,		// i
    		const int*		n,		// i
			const double*	alpha,	// i
			const double*	x,		// i[]
			const int*		incx,	// i
			const double*	y,		// i[]
			const int*		incy,	// i
			double*			A,		// io[]
			const int*		lda		// i
    );

    // y <- alpha A(') x + beta y
    void sgemv_(
    		const char*		trans,	// i 'n' or 't'
    		const int*		m,		// i
    		const int*		n,		// i
    		const float*	alpha,	// i
    		const float*	A,		// i[]
    		const int*		lda,	// i
    		const float*	x,		// i[]
    		const int*		incx,	// i
			const float*	beta,	// i
			float*			y,		// io[]
			const int*		incy	// i
    );
    void dgemv_(
    		const char*		trans,	// i 'n' or 't'
    		const int*		m,		// i
    		const int*		n,		// i
    		const double*	alpha,	// i
    		const double*	A,		// i[]
    		const int*		lda,	// i
    		const double*	x,		// i[]
    		const int*		incx,	// i
			const double*	beta,	// i
			double*			y,		// io[]
			const int*		incy	// i
    );

    // C <- alpha A(') B(') + beta C
    void dgemm_(
    		const char*		transA,	// i	'n' or 't' of A
    		const char* 	transB,	// i	'n' or 't' of B
    		const int*		m,		// i	C:mxn
    		const int*		n,		// i	C:mxn
    		const int*		k,		// i	A:mxk, B:kxn
    		const double*	alpha,	// i
    		const double*	A,		// i[]
    		const int*		lda,	// i
    		const double*	B,		// i[]
    		const int*		ldb,	// i
    		const double*	beta,	// i
    		double*			C,		// io[]
    		const int*		ldc		// i
    );
    void sgemm_(
    		const char*		transA,	// i	'n' or 't' of A
    		const char* 	transB,	// i	'n' or 't' of B
    		const int*		m,		// i	C:mxn
    		const int*		n,		// i	C:mxn
    		const int*		k,		// i	A:mxk, B:kxn
    		const float*	alpha,	// i
    		const float*	A,		// i[]
    		const int*		lda,	// i
    		const float*	B,		// i[]
    		const int*		ldb,	// i
    		const float*	beta,	// i
    		float*			C,		// io[]
    		const int*		ldc		// i
    );


    /* Lapack part*/
	// solve A*x = B
    void dgesv_(
    		const int*		n,		// i
    		const int*		nrhs,	// i
    		double*			A,		// io[]	A = P*L*U
    		const int*		lda,	// i
    		int*			ipiv,	// o[]	permutation matrix P
    		double*			B,		// io[]
    		const int*		ldb,	// i
    		int*			info	// o	= 0: successful exit
									//		< 0: if info == -i, the i-th argument had an illegal value
									//		> 0: if info == i,  U(i,i) is exactly zero.
    );
    void sgesv_(
    		const int*		n,		// i
    		const int*		nrhs,	// i
    		float*			A,		// io[]	A = P*L*U
    		const int*		lda,	// i
    		int*			ipiv,	// o[]	permutation matrix P
    		float*			B,		// io[]
    		const int*		ldb,	// i
    		int*			info	// o	= 0: successful exit
									//		< 0: if info == -i, the i-th argument had an illegal value
									//		> 0: if info == i,  U(i,i) is exactly zero.
    );
    // solve A*x = B where A is positive definite and diagonally dominant
    void dposv_(
    		const char*		uplo,	// i	'u' or 'l'
    		const int*		n,		// i	number of linear equations
    		const int*		nrhs,	// i	number of rhs columns
    		double*			A,		// io[]
    		const int*		lda,	// i
    		const double*	B,		// io[]	rhs and solution if info == 0
    		const int*		ldb,	// i
    		int*			info	// o	= 0: successful exit
    								//		< 0: if info == -i, the i-th argument had an illegal value
    								//		> 0: if info == i,  the leading minor of order i is
    								//				not positive definite, the factorization is not done.
    );
    void sposv_(
    		const char*		uplo,	// i	'u' or 'l'
    		const int*		n,		// i	number of linear equations
    		const int*		nrhs,	// i	number of rhs columns
    		float*			A,		// io[]
    		const int*		lda,	// i
    		const float*	B,		// io[]	rhs and solution if info == 0
    		const int*		ldb,	// i
    		int*			info	// o	= 0: successful exit
    								//		< 0: if info == -i, the i-th argument had an illegal value
    								//		> 0: if info == i,  the leading minor of order i is
    								//				not positive definite, the factorization is not done.
    );

    // solve A*x = B where A is symmetric but not necessary positive definite
    void dsysv_(
    		const char*		uplo,	// i	'u' or 'l'
    		const int*		n,		// i	number of linear equations
    		const int*		nrhs,	// i	number of rhs columns
    		const double*	A,		// io[] (lda, n)
    		const int*		lda,	// i
    		int*			ipiv,	// o[]
    		double*			B,		// io[]	(ldb, nrhs) rhs and solution if info == 0
    		const int*		ldb,	// i[]
    		double*			work,	// t[]	work[0] shows the optimal size when info == 0 on exit
    		int*			lwork,	// i	lwork>=n*NB, the optimal blocksize of dsytrf
    		int*			info	// o	= 0: successful exit
    		  	  	  	  	  	  	//		< 0: if INFO = -i, the i-th argument had an illegal value
    		  	  	  	  	  	    //		> 0: if INFO = i, D(i,i) is exactly zero.  The factorization has
    		  	  	  	  	  	  	//			 been completed, but the block	diagonal matrix	D is exactly singu-
    		  	  	  	  	  	  	//			 lar, so the solution could not be computed.
    );


    void ssysv_(
    		const char*		uplo,	// i	'u' or 'l'
    		const int*		n,		// i	number of linear equations
    		const int*		nrhs,	// i	number of rhs columns
    		const float*	A,		// io[] (lda, n)
    		const int*		lda,	// i
    		int*			ipiv,	// o[]
    		float*			B,		// io[]	(ldb, nrhs) rhs and solution if info == 0
    		const int*		ldb,	// i[]
    		float*			work,	// t[]	work[0] shows the optimal size when info == 0 on exit
    		int*			lwork,	// i	lwork>=n*NB, the optimal blocksize of dsytrf
    		int*			info	// o	= 0: successful exit
    		  	  	  	  	  	  	//		< 0: if INFO = -i, the i-th argument had an illegal value
    		  	  	  	  	  	    //		> 0: if INFO = i, D(i,i) is exactly zero.  The factorization has
    		  	  	  	  	  	  	//			 been completed, but the block	diagonal matrix	D is exactly singu-
    		  	  	  	  	  	  	//			 lar, so the solution could not be computed.
    );

}

// wrapper function
// performance test: not done
// easier to use due to the constant reference
// can be used as sdot(2, a, 1, b, 1) from outside
inline float sdot(
		const int& 		n, 		// dimension
		const float* 	dx, 	// []vector x
		const int& 		incx, 	// index increment of each access of x
		const float* 	dy, 	// []vector y
		const int& 		incy	// index increment of each access of y
){
	return sdot_(&n, dx, &incx, dy, &incy);
}

// wrapper function of sgesv_
inline void sgesv(
		const int&		n,		// i
		const int&		nrhs,	// i
		float*			A,		// io[]	A = P*L*U
		const int&		lda,	// i
		int*			ipiv,	// o[]	permutation matrix P
		float*			B,		// io[]
		const int&		ldb,	// i
		int&			info	// o	= 0: successful exit
								//		< 0: if info == -i, the i-th argument had an illegal value
								//		> 0: if info == i,  U(i,i) is exactly zero.
){
	sgesv_(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
}

} // end of namespace caplet


#endif /* CAPLET_BLAS_H_ */
