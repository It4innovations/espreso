#ifndef SOLID45EQNUMBGLOBAL_H_
#define SOLID45EQNUMBGLOBAL_H_
#include "utility.h"
#include "StiffnessLocal.h"

class CSolid45EqNumbGlobal {
public:
	longInt ieq[24];
	CStiffnessLocal * stif_loc;
public:
	CSolid45EqNumbGlobal();
	virtual ~CSolid45EqNumbGlobal();
};

#endif /* SOLID45EQNUMBGLOBAL_H_ */
