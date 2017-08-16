
#include "instance.h"

#include "../mesh/structures/mesh.h"
#include "../solver/generic/SparseMatrix.h"
#include "solution.h"

using namespace espreso;

Instance::Instance(const Mesh &mesh)
: domains(mesh.parts()),
  domainDOFCount(_domainDOFCount),
  properties(_properties),
  neighbours(mesh.neighbours()),
  clustersMap(mesh.getContinuityPartition()),
  K(_K), N1(_N1), N2(_N2), RegMat(_RegMat),
  M(_M),
  R(_R), f(_f),
  B0(_B0),
  B0subdomainsMap(_B0subdomainsMap),
  B1(_B1),
  B1subdomainsMap(_B1subdomainsMap),
  B1clustersMap(_B1clustersMap),
  B1c(_B1c), LB(_LB), B1duplicity(_B1duplicity),
  inequality(_inequality), inequalityC(_inequalityC),
  block(_block)
{
	domainDOFCount.resize(domains);

	K.resize(domains);
	N1.resize(domains);
	N2.resize(domains);
	RegMat.resize(domains);

	M.resize(domains);
	R.resize(domains);
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

	computeKernelCallback = [] (REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernel is empty function. Fill it in assembler.";
	};

	computeKernelsCallback = [] (REGULARIZATION regularization, size_t scSize, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernels is empty function. Fill it in assembler.";
	};

	assembleB0Callback = [] (B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: assembleB0 is empty function. Fill it in assembler.";
	};
}

Instance::Instance(Instance &other, Matrices &share)
: domains(other.domains),
  domainDOFCount(share & Matrices::K ? other.domainDOFCount :_domainDOFCount),
  properties(other.properties),
  neighbours(other.neighbours),
  clustersMap(other.clustersMap),

  // shared K -> also share kernels and regularization matrix
  K(share & Matrices::K ? other.K : _K),
  N1(share & Matrices::K ? other.N1 : _N1),
  N2(share & Matrices::K ? other.N2 : _N2),
  RegMat(share & Matrices::K ? other.RegMat : _RegMat),

  M(share & Matrices::M ? other.M : _M),
  R(share & Matrices::R ? other.R : _R),
  f(share & Matrices::f ? other.f : _f),

  B0(share & Matrices::B0 ? other.B0 : _B0),
  B0subdomainsMap(share & Matrices::B0 ? other.B0subdomainsMap : _B0subdomainsMap),

  B1(share & Matrices::B1 ? other.B1 : _B1),
  B1subdomainsMap(share & Matrices::B1 ? other.B1subdomainsMap : _B1subdomainsMap),
  B1clustersMap(share & Matrices::B1 ? other.B1clustersMap : _B1clustersMap),
  B1c(share & Matrices::B1c ? other.B1c : _B1c),
  LB(share & Matrices::B1 ? other.LB : _LB),
  B1duplicity(share & (Matrices::B1 | Matrices::K) ? other.B1duplicity : _B1duplicity),
  inequality(share & Matrices::B1 ? other.inequality : _inequality),
  inequalityC(share & Matrices::B1 ? other.inequalityC : _inequalityC),
  block(share & Matrices::B1 ? other.block : _block)
{
	if (!(share & Matrices::K)) {
		domainDOFCount.resize(domains);
		K.resize(domains);
		N1.resize(domains);
		N2.resize(domains);
		RegMat.resize(domains);
	}

	if (!(share & Matrices::M)) {
		M.resize(domains);
	}
	if (!(share & Matrices::R)) {
		R.resize(domains);
	}
	if (!(share & Matrices::f)) {
		f.resize(domains);
	}

	if (!(share & Matrices::B0)) {
		B0.resize(domains);
		B0subdomainsMap.resize(domains);

		for (size_t d = 0; d < domains; d++) {
			B0[d].rows = 0;
			B0[d].cols = 0;
			B0[d].nnz = 0;
			B0[d].type = 'G';
		}
	}

	if (!(share & Matrices::B1)) {
		B1.resize(domains);
		B1subdomainsMap.resize(domains);
		LB.resize(domains);

		inequality.resize(domains);
		inequalityC.resize(domains);

		block.resize(3);

		for (size_t d = 0; d < domains; d++) {
			B1[d].rows = 0;
			B1[d].cols = 0;
			B1[d].nnz = 0;
			B1[d].type = 'G';
		}
	}

	if (!(share & Matrices::B1c)) {
		B1c.resize(domains);
	}

	if (!((share & Matrices::B1c) & (share & Matrices::K))) {
		B1duplicity.resize(domains);
	}

	computeKernelCallback = [] (REGULARIZATION regularization, size_t scSize, size_t domain, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernel is empty function. Fill it in assembler.";
	};

	computeKernelsCallback = [] (REGULARIZATION regularization, size_t scSize, bool ortogonalCluster) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: computeKernels is empty function. Fill it in assembler.";
	};

	assembleB0Callback = [] (B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: assembleB0 is empty function. Fill it in assembler.";
	};
}

void Instance::clear()
{
	K.clear();
	K.resize(domains);
	N1.clear();
	N1.resize(domains);
	N2.clear();
	N2.resize(domains);
	RegMat.clear();
	RegMat.resize(domains);

	M.clear();
	M.resize(domains);
	R.clear();
	R.resize(domains);
	f.clear();
	f.resize(domains);

	B0.clear();
	B0.resize(domains);
	B0subdomainsMap.clear();
	B0subdomainsMap.resize(domains);

	B1.clear();
	B1.resize(domains);
	B1subdomainsMap.clear();
	B1subdomainsMap.resize(domains);
	B1duplicity.clear();
	B1duplicity.resize(domains);
	B1c.clear();
	B1c.resize(domains);
	LB.clear();
	LB.resize(domains);

	inequality.clear();
	inequality.resize(domains);
	inequalityC.clear();
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



