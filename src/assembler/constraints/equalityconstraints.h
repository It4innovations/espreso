
#ifndef SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include <vector>
#include <cstddef>

namespace espreso {

class Element;
class Region;
class Mesh;
enum class Property;
class Instance;
struct Step;
class SparseMatrix;

struct EqualityConstraints
{
	EqualityConstraints(
			Instance &instance,
			const Mesh &mesh,
			const std::vector<Element*> &gluedElements,
			const std::vector<Element*> &gluedInterfaceElements,
			const std::vector<Property> &gluedDOFs,
			const std::vector<size_t> &gluedDOFsMeshOffsets,
			bool interfaceElementContainsGluedDOFs = false);

	void insertDirichletToB1(const Step &step, bool withRedundantMultiplier);
	void insertElementGluingToB1(const Step &step, bool withRedundantMultiplier, bool withScaling);

	void insertCornersGluingToB0();
	void insertKernelsGluingToB0(const std::vector<SparseMatrix> &kernels);

protected:
	std::vector<esglobal> computeLambdasID(const Step &step, bool withRedundantMultiplier);

	Instance &_instance;
	const Mesh &_mesh;
	const std::vector<Element*> &_gluedElements;
	const std::vector<Element*> &_gluedInterfaceElements;
	const std::vector<Property> &_gluedDOFs;
	const std::vector<size_t> &_gluedDOFsMeshOffsets;
	bool _interfaceElementContainsGluedDOFs;
};

}



#endif /* SRC_ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
