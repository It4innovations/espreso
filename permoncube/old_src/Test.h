/*
 * Test.h
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#ifndef TEST_H_
#define TEST_H_

#include "utility.h"

class CTest {
public:
	CTest();
	virtual ~CTest();

public:
	static void store1darray2ascii(int *i, double * x, int n);
};

#endif /* TEST_H_ */
