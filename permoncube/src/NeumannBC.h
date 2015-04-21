/*
 * NeumannBC.h
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#ifndef NEUMANNBC_H_
#define NEUMANNBC_H_

#include "utility.h"

class CNeumannBC {
public:
	  int *ind;
	  double *val;
	  int n;

public:
	CNeumannBC();
	virtual ~CNeumannBC();
};

#endif /* NEUMANNBC_H_ */
