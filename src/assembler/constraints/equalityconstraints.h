
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include "constraints.h"

namespace espreso {

class Element;
class Region;
enum class Property;
class Instance;
struct Step;

struct EqualityConstraints
{
	static void insertDirichletToB1(Constraints &constraints, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs);
	static void insertElementGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &K);
	static void insertMortarGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);

	static void insertDomainGluingToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	static void insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &kernel);
	static void insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Element*> &DOFs, const std::vector<SparseMatrix> &kernel);

	static void insertDirichletToB1(Instance &instance, const Step &step, const std::vector<Region*> &regions, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs, const std::vector<size_t> &DOFsOffsets, bool withRedundantMultiplier);
	static void insertElementGluingToB1(Instance &instance, const Step &step, const std::vector<int> &neighbours, const std::vector<Region*> &regions, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<size_t> &DOFsOffsets, bool withRedundantMultiplier, bool withScaling);

	static void insertCornersGluingToB0(Instance &instance, const std::vector<Element*> &elements, const std::vector<size_t> &DOFsOffsets);
	static void insertKernelsGluingToB0(Instance &instance, const std::vector<Element*> &elements, const std::vector<Element*> &DOFsSource, const std::vector<size_t> &DOFsOffsets, const std::vector<SparseMatrix> &kernels, bool getDOFsSourceIndicesFromElement = false);

protected:
	static std::vector<esglobal> computeLambdasID(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs);
	static std::vector<esglobal> computeLambdasID(Instance &instance, const Step &step, const std::vector<int> &neighbours, const std::vector<Region*> &regions, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<size_t> &DOFsOffsets, bool withRedundantMultiplier);

};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
