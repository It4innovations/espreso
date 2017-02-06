
#include "instance.h"

#include "../solver/generic/SparseMatrix.h"
#include "solution.h"

using namespace espreso;

Instance::Instance(size_t domains): domains(domains)
{
	DOFs.resize(domains);

	K.resize(domains);
	R1.resize(domains);
	R2.resize(domains);
	RegMat.resize(domains);
	f.resize(domains);

	B0.resize(domains);
	B0subdomainsMap.resize(domains);

	B1.resize(domains);
	B1subdomainsMap.resize(domains);
	B1duplicity.resize(domains);
	B1c.resize(domains);
	LB.resize(domains);

	inequality.resize(domains);
	inequalityC.resize(domains);

	block.resize(3);

	for (size_t d = 0; d < domains; d++) {
		B0[d].rows = 0;
		B0[d].cols = 0;
		B0[d].nnz = 0;
		B0[d].type = 'G';

		B1[d].rows = 0;
		B1[d].cols = 0;
		B1[d].nnz = 0;
		B1[d].type = 'G';
	}
}

Instance::~Instance()
{
	for (size_t i = 0; i < solutions.size(); i++) {
		delete solutions[i];
	}
}



