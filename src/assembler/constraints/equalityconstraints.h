
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include "constraints.h"

namespace espreso {

struct EqualityConstraints
{
	static void insertDirichletToB1(Constraints &constraints, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs);
	static void insertElementGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &K);
	static void insertMortarGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

	static void insertDomainGluingToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	static void insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &kernel);
	static void insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Element*> &DOFs, const std::vector<SparseMatrix> &kernel);

protected:
	static std::vector<esglobal> computeLambdasID(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
