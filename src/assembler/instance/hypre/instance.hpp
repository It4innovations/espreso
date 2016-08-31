
#include "instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void HypreInstance<TConstrains, TPhysics>::init()
{

	std::cout << "FILL MESS\n";
	const std::vector<Element*> &elements = _mesh.elements();

	for (size_t i = 0; i < elements.size(); i++) {
		std::cout << "cluster index: ";
		for (size_t j = 0; j < elements[i]->nodes(); j++) {
			std::cout << elements[i]->node(j) << " ";
		}
		std::cout << "\n";
		std::cout << "global index: ";
		for (size_t j = 0; j < elements[i]->nodes(); j++) {
			std::cout << _mesh.coordinates().globalIndex(elements[i]->node(j)) << " ";
		}
		std::cout << "\n";
	}

	std::cout << "FILL BC\n";

	const std::vector<Element*> &nodes = _mesh.nodes();

	for (size_t i = 0; i < nodes.size(); i++) {
		if (nodes[i]->settings().isSet(Property::TEMPERATURE)) {
			const Point &p = _mesh.coordinates()[i];
			std::cout << "node " << i << " has fixed temperature " << nodes[i]->settings(Property::TEMPERATURE).back()->evaluate(0) << "\n";
		}
	}

	 DenseMatrix Ke;
	 std::vector<double> fe;

	 for (size_t i = 0; i < elements.size(); i++) {
		 _physics.assembleStiffnessMatrix(elements[i], Ke, fe);
	 }
}

template <class TConstrains, class TPhysics>
void HypreInstance<TConstrains, TPhysics>::solve(std::vector<std::vector<double> > &solution)
{
	// SOLVER HYPRE
	int status;
	feiPtr.solve(&status);


	std::cout << "SAVE SOLUTION\n";

	solution.resize(_mesh.parts());
	for (size_t p = 0; p < _mesh.parts(); p++) {
		solution[p].resize(_mesh.coordinates().localSize(p));
	}
}

}
