
#ifndef ESPRESO_WRAPPER_H_
#define ESPRESO_WRAPPER_H_

#include <iostream>
#include <list>

#include "esconfig.h"
#include "esassemblers.h"
#include "esinput.h"
#include "esbasis.h"

struct FETI4IStructMatrix {
	FETI4IStructMatrix(eslocal offset): offset(offset), K(0, 0) { };

	std::vector<std::vector<eslocal> > eIndices;
	std::vector<std::vector<double> > eMatrices;

	// used in case of 1 subdomain
	espreso::SparseVVPMatrix<eslocal> K;

	eslocal offset;
};

struct FETI4IStructInstance {
	FETI4IStructInstance(espreso::LinearElasticity<espreso::API> data, espreso::Mesh *mesh): data(data), mesh(mesh), K(0, 0) { };
	~FETI4IStructInstance() { delete mesh; }

	espreso::LinearElasticity<espreso::API> data;
	espreso::Mesh *mesh;
	espreso::SparseCSRMatrix<eslocal> K;
};

namespace espreso {

struct DataHolder {
	static std::list<FETI4IStructMatrix*> matrices;
	static std::list<FETI4IStructInstance*> instances;
	static TimeEval timeStatistics;
};

}


#endif /* ESPRESO_WRAPPER_H_ */
