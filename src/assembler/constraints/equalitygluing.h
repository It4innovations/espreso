
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYGLUING_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYGLUING_H_

#include "constraints.h"

namespace espreso {

class EqualityGluing: public Constraints
{
public:
	EqualityGluing(Mesh &mesh): Constraints(mesh) {};

	void insertDirichletToB1(const std::vector<Element*> &nodes, const Coordinates &coordinates, const std::vector<Property> &DOFs);
	void insertElementGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	void insertMortarGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

	void insertDomainGluingToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

protected:
	std::vector<esglobal> computeLambdasID(const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYGLUING_H_ */
