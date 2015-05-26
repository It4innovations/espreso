/*
 * Test.cpp
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#include "Test.h"

CTest::CTest() {
	// TODO Auto-generated constructor stub

}

CTest::~CTest() {
	// TODO Auto-generated destructor stub
}

void CTest::store1darray2ascii(int *ind,double *x, int n){

	char filename[128];
	sprintf(filename, "%s","data/file1.txt");
	ofstream f_001(filename);
	for (int i = 0; i < n; i++) {
		f_001 <<ind[i] <<" "<<x[i] << endl;
	}
	f_001.close();

}

