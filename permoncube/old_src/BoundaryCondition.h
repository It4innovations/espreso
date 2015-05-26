
#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include "utility.h"
#include "DirichletBC.h"
#include "NeumannBC.h"
#include "ContactBC.h"

class CBoundaryCondition {
public:
	  CDirichletBC *dirBCSub;
	  CDirichletBC *dirBC_global;
	  CNeumannBC *neuBC;
	  CContactBC *conBCSub;

public:
	CBoundaryCondition();
  CBoundaryCondition (const CBoundaryCondition & _CBoundaryCondition);

	virtual ~CBoundaryCondition();
};

#endif /* BOUNDARYCONDITION_H_ */
