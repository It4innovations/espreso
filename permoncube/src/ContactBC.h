
#ifndef CONTACTBC_H_
#define CONTACTBC_H_

#include "utility.h"

class CContactBC{
public:
	  int *ind;
	  double *val;
	  int n;

public:
	  CContactBC();
	virtual ~CContactBC();
};

#endif /* CONTACTBC_H_ */
