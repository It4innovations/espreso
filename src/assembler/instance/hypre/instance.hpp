
#include "instance.h"

namespace espreso {

template <class TConstrains, class TPhysics>
void HypreInstance<TConstrains, TPhysics>::init()
{
	_physics.prepareMeshStructures();

	if (config::output::SAVE_PROPERTIES) {
		_physics.saveMeshProperties(_store);
	}

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

	std::vector<Property> DOFs = _physics.pointDOFs;

	std::cout << DOFs;

	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (nodes[i]->settings().isSet(DOFs[dof])) {
				const Point &p = _mesh.coordinates()[i];
				std::cout << "node " << i << " has fixed " << DOFs[dof] << " to " << nodes[i]->settings(DOFs[dof]).back()->evaluate(i) << "\n";
			}
		}
	}

	 DenseMatrix Ke;
	 std::vector<double> fe;

	 const std::vector<eslocal> &partition = _mesh.getPartition();

	 for (size_t subdomain = 0; subdomain < _mesh.parts(); subdomain++) {
		 for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

			 // compute element matrix
			 _physics.assembleStiffnessMatrix(elements[e], Ke, fe);

			 std::cout << "Element matrix: \n";
			 for (size_t nx = 0; nx < elements[e]->nodes(); nx++) {
				for (size_t dx = 0; dx < DOFs.size(); dx++) {
					size_t row = nodes[elements[e]->node(nx)]->DOFIndex(subdomain, dx);
					std::cout << "row:";
					for (size_t ny = 0; ny < elements[e]->nodes(); ny++) {
						for (size_t dy = 0; dy < DOFs.size(); dy++) {
							size_t column = nodes[elements[e]->node(ny)]->DOFIndex(subdomain, dy);
							std::cout << " " << std::setprecision(1) << Ke(dx * elements[e]->nodes() + nx, dy * elements[e]->nodes() + ny);
						}
					}
					std::cout << " rhs: "<< fe[dx * elements[e]->nodes() + nx] << "\n";
				}
			}

		 }
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

	if (config::output::SAVE_RESULTS) {
		_physics.saveMeshResults(_store, solution);
	}
}

}
