/*
 * LinearAlgebra.h
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#ifndef LINEARALGEBRA_H_
#define LINEARALGEBRA_H_

#include "utility.h"
#include "Solid45EqNumbGlobal.h"
#include "Fem.h"
#include "ElementLibrary.h"
#include "KSparse.h"

class CSolid45EqNumbGlobal;
class CFem;
class CKSparce;

class CLinearAlgebra {
private:
	CLinearAlgebra();
	virtual ~CLinearAlgebra();

public:

	static double inverse_matrix_3x3(double *J,double *iJ);

	static int compare (const void * a, const void * b);
	static int compare2 (const void * a, const void * b);
	static int compareDouble(const void * a, const void * b);
	static int compareInt(const void * a, const void * b);
	static int compareLongInt(const void * a, const void * b);
  static bool compare_couple(const void * a, const void * b); 
  static bool equal_couple(const void * a, const void * b); 

	static double norm_v(double * x,const int &n);
	static void add_vec_a2b(double *a,double *b,double ka,double kb,const int &n);
//	static void zeroing(double *x ,const int &n);



};

#endif /* LINEARALGEBRA_H_ */
