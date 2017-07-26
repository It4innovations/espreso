
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>
#include <list>
#include <vector>

#include "../../basis/logging/timeeval.h"

namespace espreso {
struct Environment;
class Instance;
class Physics;
class SolverBase;
class APIMesh;
class OutputConfiguration;
class ESPRESOSolver;
class LinearSolver;

namespace output {
class ResultStoreList;
}
}

struct FETI4IStructMatrix {
	FETI4IStructMatrix(eslocal type, eslocal offset): type(type), offset(offset) {};

	std::vector<eslocal> eType;
	std::vector<std::vector<eslocal> > eNodes;
	std::vector<std::vector<eslocal> > eDOFs;
	std::vector<std::vector<double> > eMatrices;

	eslocal type;
	eslocal offset;
};

struct FETI4IStructInstance {
	FETI4IStructInstance(FETI4IStructMatrix &matrix, eslocal *l2g, size_t size);
	~FETI4IStructInstance();

	espreso::Instance *instance;
	espreso::Physics * physics;
	espreso::LinearSolver *linearSolver;
	espreso::output::ResultStoreList *store;
	espreso::SolverBase *solver;

	espreso::APIMesh *mesh;
	espreso::OutputConfiguration *output;
	espreso::ESPRESOSolver *configuration;
};

namespace espreso {

struct DataHolder {
	static Environment environment;
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructInstance*> instances;
	static TimeEval timeStatistics;
};

}


#endif /* ESPRESO_WRAPPER_H_ */
