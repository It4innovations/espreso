
#ifndef SRC_APP_APITEST_H_
#define SRC_APP_APITEST_H_

#include <config/ecf/ecf.h>

#include "factory/factory.h"

#include "basis/matrices/denseMatrix.h"
#include "mesh/mesh.h"
#include "assembler/step.h"
#include "assembler/instance.h"
#include "assembler/constraints/constraints.h"
#include "assembler/physics/physics.h"
#include "assembler/physicssolver/assembler.h"
#include "mesh/store/elementstore.h"
#include "output/result/resultstore.h"
#include "solver/generic/SparseMatrix.h"

namespace espreso {

class APITestESPRESODataProvider {

public:
	APITestESPRESODataProvider(int *argc, char ***argv): ecf(argc, argv), mesh(ecf), factory(ecf, mesh, store)
	{
		if (factory._loadSteps.size() > 1) {
			ESINFO(GLOBAL_ERROR) << "APITEST: Cannot test instance with more loadsteps.";
		}
		factory._loader->_physics[0]->preprocessData();
	}

	int matrixType()
	{
		return static_cast<int>(factory._loader->_assemblers[0]->physics.getMatrixType(0));
	}

	size_t elements()
	{
		return factory._mesh->elements->size;
	}

	void addElementMatrix(FETI4IMatrix &K, size_t e)
	{
//		FETI4IInt eType = static_cast<int>(factory._mesh->_elems->epointers->datatarray()[e]->type);
//		std::vector<FETI4IInt> nodes;
//		for (size_t n = 0; n < factory._mesh->elements()[e]->nodes(); n++) {
//			nodes.push_back(factory._mesh->elements()[e]->node(n));
//		}
//		std::vector<FETI4IInt> dofs;
//		factory._loader->_physics[0]->fillDOFsIndices(factory._mesh->elements()[e], 0, dofs);
//		DenseMatrix eK, eM, eR, eF;
//		factory._loader->_physics[0]->updateMatrix(step, Matrices::K, factory._mesh->elements()[e], eK, eM, eR, eF, factory._loader->_instances[0]->solutions);
//		FETI4IAddElement(K, eType, nodes.size(), nodes.data(), dofs.size(), dofs.data(), eK.values());
	}

	void computeRHS(std::vector<FETI4IReal> &rhs)
	{
//		rhs.resize(factory._loader->_instances[0]->domainDOFCount[0]);
//		std::vector<FETI4IInt> dofs;
//		for (size_t e = 0; e < elements(); e++) {
//			factory._loader->_physics[0]->fillDOFsIndices(factory._mesh->elements()[e], 0, dofs);
//			DenseMatrix eK, eM, eR, eF;
//			factory._loader->_physics[0]->updateMatrix(step, Matrices::f, factory._mesh->elements()[e], eK, eM, eR, eF, factory._loader->_instances[0]->solutions);
//			for (size_t dof = 0; dof < dofs.size(); dof++) {
//				rhs[dofs[dof]] = eF(0, dof);
//			}
//		}
	}

	void fillL2G(std::vector<FETI4IInt> &l2g)
	{
//		for (size_t n = 0; n < factory._mesh->nodes().size(); n++) {
//			if (factory._mesh->nodes()[n]->parentElements().size()) {
//				l2g.push_back(factory._mesh->coordinates().globalIndex(n));
//			}
//		}
	}

	void fillDirichlet(std::vector<FETI4IInt> &dirichlet_indices, std::vector<FETI4IReal> &dirichlet_values)
	{
		factory._loader->_physics[0]->_constraints->B1DirichletInsert(step);
		for (auto it = factory._loader->_instances[0]->B1[0].J_col_indices.begin(); it != factory._loader->_instances[0]->B1[0].J_col_indices.end(); ++it) {
			dirichlet_indices.push_back(*it - 1);
		}
		dirichlet_values.insert(dirichlet_values.end(), factory._loader->_instances[0]->B1c[0].begin(), factory._loader->_instances[0]->B1c[0].end());
	}

	void fillNeighbours(std::vector<FETI4IMPIInt> &neighbors)
	{
		neighbors = factory._mesh->neighbors;
	}

protected:
	ECF ecf;
	Mesh mesh;
	ResultStore store;
	Factory factory;
	Step step;
};

}



#endif /* SRC_APP_APITEST_H_ */
