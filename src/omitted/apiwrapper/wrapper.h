
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <config/ecf/ecf.h>
#include <iostream>
#include <list>
#include <vector>

#include "basis/logging/timeeval.h"

namespace espreso {
struct Environment;
class DataHolder;
class Physics;
class SubStepSolver;
class LoadStepSolver;
class Assembler;
class Mesh;
class ECF;
class FETISystemSolver;
class ResultStore;
}

struct FETI4IStructMatrix {
	FETI4IStructMatrix(esint type, esint offset): type(type), offset(offset) {};

	std::vector<esint> eType;
	std::vector<std::vector<esint> > eNodes;
	std::vector<std::vector<esint> > eDOFs;
	std::vector<std::vector<double> > eMatrices;

	esint type;
	esint offset;
};

struct FETI4IStructInstance {
	FETI4IStructInstance(FETI4IStructMatrix &matrix, esint *l2g, size_t size);
	~FETI4IStructInstance();

	espreso::DataHolder *instance;
	espreso::Physics * physics;
	espreso::FETISystemSolver *linearSolver;
	espreso::ResultStore *store;
	espreso::Step *step;
	espreso::Assembler *assembler;
	espreso::SubStepSolver *timeStepSolver;
	espreso::LoadStepSolver *loadStepSolver;

	espreso::Mesh *mesh;
	espreso::ECF configuration;
};

namespace espreso {

struct APIDataHolder {
	static ECF *configuration;
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructInstance*> instances;
	static TimeEval timeStatistics;
};

}


#endif /* ESPRESO_WRAPPER_H_ */
