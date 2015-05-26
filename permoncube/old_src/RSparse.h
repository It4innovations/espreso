 
#ifndef RSPARSE_H_
#define RSPARSE_H_

#include "utility.h"
#include "Fem.h"
#include "mpi.h"

class CFem;

class CRSparse {

public:
  double *Rfull;
  int n_row;
  int n_col;
	CRSparse(CFem *fem, int neqSub);
	virtual ~CRSparse();
};

#endif /* BSPARSE_H_ */ 
 
