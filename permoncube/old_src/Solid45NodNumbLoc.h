/*
 * Solid45NodNumbLoc.h
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#ifndef SOLID45NODNUMBLOC_H_
#define SOLID45NODNUMBLOC_H_

#include "utility.h"

class CSolid45NodNumbLoc {
public:
	int inod_loc[8];
public:
	CSolid45NodNumbLoc();
	virtual ~CSolid45NodNumbLoc();
};

#endif /* SOLID45NODNUMBLOC_H_ */
