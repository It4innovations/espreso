
#ifndef DIRICHLETBC_H_
#define DIRICHLETBC_H_

#include "utility.h"

class CDirichletBC {
public:
	  int *ind;
	  double *val;
	  int n;

public:
	  CDirichletBC();
	virtual ~CDirichletBC();
};

#endif /* DIRICHLETBC_H_ */
