
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include "constraints.h"

namespace espreso {

class Element;
class Region;
enum class Property;
class NewInstance;

struct EqualityConstraints
{
	static void insertDirichletToB1(Constraints &constraints, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs);
	static void insertElementGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &K);
	static void insertMortarGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

	static void insertDomainGluingToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	static void insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &kernel);
	static void insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Element*> &DOFs, const std::vector<SparseMatrix> &kernel);

	static void insertDirichletToB1(NewInstance &instance, const std::vector<Region*> &regions, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs, bool withRedundantMultiplier);
	static void insertElementGluingToB1(NewInstance &instance, const std::vector<int> &neighbours, const std::vector<Region*> &regions, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, bool withRedundantMultiplier, bool withScaling);

	static void insertCornersGluingToB0(NewInstance &instance, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	static void insertKernelsGluingToB0(NewInstance &instance, const std::vector<Element*> &elements, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs);

protected:
	static std::vector<esglobal> computeLambdasID(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	static std::vector<esglobal> computeLambdasID(NewInstance &instance, const std::vector<int> &neighbours, const std::vector<Region*> &regions, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, bool withRedundantMultiplier);

};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
