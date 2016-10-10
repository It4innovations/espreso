
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include "constraints.h"

namespace espreso {

class EqualityConstraints: public Constraints
{
public:
	EqualityConstraints(Mesh &mesh): Constraints(mesh) {};

	void insertDirichletToB1(const std::vector<Element*> &nodes, const std::vector<Property> &DOFs);
	void insertElementGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	void insertMortarGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

	void insertDomainGluingToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	void insertKernelsToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &kernel);
	void insertKernelsToB0(const std::vector<Element*> &elements, const std::vector<Element*> &DOFs, const std::vector<SparseMatrix> &kernel);

protected:
	std::vector<esglobal> computeLambdasID(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
