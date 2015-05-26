/*
 * LinearAlgebra.cpp
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#include "LinearAlgebra.h"

CLinearAlgebra::CLinearAlgebra() {
	// TODO Auto-generated constructor stub

}

CLinearAlgebra::~CLinearAlgebra() {
	// TODO Auto-generated destructor stub
}

double CLinearAlgebra::inverse_matrix_3x3(double *A, double *iA) {

	double determinant = +A[0] * (A[4] * A[8] - A[5] * A[7])
			- A[3] * (A[1] * A[8] - A[7] * A[2])
			+ A[6] * (A[1] * A[5] - A[4] * A[2]);
	double invdet = 1. / determinant;

	iA[0] = (A[4] * A[8] - A[5] * A[7]) * invdet;
	iA[1] = -(A[3] * A[8] - A[6] * A[5]) * invdet;
	iA[2] = (A[3] * A[7] - A[6] * A[4]) * invdet;
	iA[3] = -(A[1] * A[8] - A[7] * A[2]) * invdet;
	iA[4] = (A[0] * A[8] - A[6] * A[2]) * invdet;
	iA[5] = -(A[0] * A[7] - A[1] * A[6]) * invdet;
	iA[6] = (A[1] * A[5] - A[2] * A[4]) * invdet;
	iA[7] = -(A[0] * A[5] - A[2] * A[3]) * invdet;
	iA[8] = (A[0] * A[4] - A[1] * A[3]) * invdet;

	return determinant;
}

bool CLinearAlgebra::compare_couple(const void * a, const void * b) {
    return ((struct CCoupleIntDouble*) a)->ind
        - ((struct CCoupleIntDouble*) b)->ind;
}

bool CLinearAlgebra::equal_couple(const void * a, const void * b) {
    return ((struct CCoupleIntDouble*) a)->ind ==
         ((struct CCoupleIntDouble*) b)->ind;
}

int CLinearAlgebra::compareInt(const void * a, const void * b) {
	return ( *(int*)a - *(int*)b );
}

int CLinearAlgebra::compareLongInt(const void * a, const void * b) {
	return ( *(longInt*)a - *(longInt*)b );
}

int CLinearAlgebra::compare(const void * a, const void * b) {
	return ((struct CCustomData*) a)->data - ((struct CCustomData*) b)->data;
}

int CLinearAlgebra::compare2(const void * a, const void * b) {
	return ((struct CCustomData*) a)->orig_pos
			- ((struct CCustomData*) b)->orig_pos;
}

int CLinearAlgebra::compareDouble(const void * a, const void * b) {
	if (*(double*) a >= *(double*) b) {
		return 1;
	} else {
		return -1;
	}
}

double CLinearAlgebra::norm_v(double *x, const int &n) {
	double sum = 0.;
	for (int i = 0; i < n; i++) {
		sum += x[i] * x[i];
	}

	return sqrt(sum);
}

void CLinearAlgebra::add_vec_a2b(double *a, double *b, double ka, double kb,
		const int &n) {
	for (int i = 0; i < n; i++) {
		b[i] = kb * b[i] + ka * a[i];
	}
}

