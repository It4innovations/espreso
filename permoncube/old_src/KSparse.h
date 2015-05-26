
#ifndef KSPARCE_H_
#define KSPARCE_H_

#include "utility.h"
#include "Fem.h"
#include "Solid45EqNumbGlobal.h"
#include "StiffnessLocal.h"
#include "BSparse.h"
#include "CoupleIntDouble.h"


class CBSparse;
class CFem;
class CCoupleIntDouble;

class CKSparse {

public:
  MPI_Comm comm;
	double *val;
	int *col_ind;
	int *row_ptr;
  int n_row;    
//  struct s_i_d{ int ind; double val; };

	CKSparse(MPI_Comm comm);
	virtual ~CKSparse();

	void multAx(double *Ax, double *x, int neqSub,int regularized);
	void fe_symbolic_form_matrix(CFem *fem, CDomain *domainG,CSolid45EqNumbGlobal *stif_glob_numbering);
	void fe_numeric_form_matrix(CSolid45EqNumbGlobal * stif_glob_numb,
															CStiffnessLocal * stif_loc,CFem * fem);
};


#endif /* KSPARCE_H_ */
