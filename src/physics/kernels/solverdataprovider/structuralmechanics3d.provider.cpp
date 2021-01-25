
#include "structuralmechanics3d.provider.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/fetidatastore.h"
#include "math/math.h"
#include "math/matrix.type.h"
#include "math/matrix.dense.feti.h"
#include "math/matrix.csr.feti.h"
#include "math/matrix.csr.distributed.h"
#include "math/vector.dense.h"
#include "math/vector.dense.distributed.h"

#include <cmath>
#include <algorithm>

#include "mkl.h"
#include "mkl_solvers_ee.h"

using namespace espreso;

MatrixType StructuralMechanics3DSolverDataProvider::General::getMatrixType()
{
	if (_configuration.type == LoadStepSolverConfiguration::TYPE::HARMONIC) {
		if (_configuration.harmonic_solver.damping.rayleigh.type != HarmonicRayleighDampingConfiguration::TYPE::NONE) {
			return MatrixType::REAL_UNSYMMETRIC;
		}
		if (_configuration.harmonic_solver.damping.coriolis_effect.coriolis_damping) {
			return MatrixType::REAL_UNSYMMETRIC;
		}
		return MatrixType::REAL_SYMMETRIC_INDEFINITE;
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void StructuralMechanics3DSolverDataProvider::General::dirichletIndices(std::vector<std::pair<esint, esint>> &indices)
{
	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		if (it->second.all.value.size() || it->second.x.value.size()) {
			for (auto i = region->nodes->datatarray().begin(); i != region->nodes->datatarray().end(); ++i) {
				indices.push_back(std::make_pair(*i, 0));
			}
		}
		if (it->second.all.value.size() || it->second.y.value.size()) {
			for (auto i = region->nodes->datatarray().begin(); i != region->nodes->datatarray().end(); ++i) {
				indices.push_back(std::make_pair(*i, 1));
			}
		}
		if (it->second.all.value.size() || it->second.z.value.size()) {
			for (auto i = region->nodes->datatarray().begin(); i != region->nodes->datatarray().end(); ++i) {
				indices.push_back(std::make_pair(*i, 2));
			}
		}
	}
}

void StructuralMechanics3DSolverDataProvider::General::dirichletValues(std::vector<double> &values)
{
	size_t offset = 0;
	double *coors = reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data());
	auto eval = [&] (Evaluator *evaluator, tarray<esint> &nodes) {
		evaluator->evalSelectedSparse(nodes.size(), nodes.data(), Evaluator::Params().coords(3, coors), values.data() + offset);
		offset += nodes.size();
	};

	auto pick = [&] (ECFExpressionOptionalVector &vector, tarray<esint> &nodes) {
		if (vector.all.value.size()) {
			eval(vector.all.evaluator, nodes);
			eval(vector.all.evaluator, nodes);
			eval(vector.all.evaluator, nodes);
		} else {
			if (vector.x.value.size()) {
				eval(vector.x.evaluator, nodes);
			}
			if (vector.y.value.size()) {
				eval(vector.y.evaluator, nodes);
			}
			if (vector.z.value.size()) {
				eval(vector.z.evaluator, nodes);
			}
		}
	};

	for (auto it = _configuration.displacement.begin(); it != _configuration.displacement.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		pick(it->second, region->nodes->datatarray());
	}
}

void StructuralMechanics3DSolverDataProvider::General::inequalityIndices(std::vector<std::pair<esint, esint>> &indices)
{
	for (auto it = _configuration.obstacle.begin(); it != _configuration.obstacle.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		for (int dof = 0; dof < 3; ++dof) {
			for (auto i = region->nodes->datatarray().begin(); i != region->nodes->datatarray().end(); ++i) {
				indices.push_back(std::make_pair(*i, dof));
			}
		}
	}
}

static void evaluate(std::vector<double> &values, std::map<std::string, ECFExpressionVector> &bc)
{
	size_t offset = 0;
	Evaluator::Params params;
	params.coords(3, reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data()));

	for (auto it = bc.begin(); it != bc.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		auto nodes = region->nodes->datatarray();

		auto eval = [&] (ECFExpression &expr) {
			expr.evaluator->evalSelectedSparse(nodes.size(), nodes.data(), params, values.data() + offset);
			offset += region->nodes->datatarray().size();
		};

		eval(it->second.x);
		eval(it->second.y);
		eval(it->second.z);
	}
}

void StructuralMechanics3DSolverDataProvider::General::inequalityNormals(std::vector<double> &values)
{
	evaluate(values, _configuration.normal_direction);
}

void StructuralMechanics3DSolverDataProvider::General::inequalityGaps(std::vector<double> &values)
{
	evaluate(values, _configuration.obstacle);
}

StructuralMechanics3DSolverDataProvider::FETI::~FETI()
{
	if (_RegMat) { delete _RegMat; }
}

MatrixType StructuralMechanics3DSolverDataProvider::FETI::getMatrixType(esint domain)
{
	if (_configuration.type == LoadStepSolverConfiguration::TYPE::HARMONIC) {
		if (_configuration.harmonic_solver.damping.rayleigh.type != HarmonicRayleighDampingConfiguration::TYPE::NONE) {
			return MatrixType::REAL_UNSYMMETRIC;
		}
		if (_configuration.harmonic_solver.damping.coriolis_effect.coriolis_damping) {
			return MatrixType::REAL_UNSYMMETRIC;
		}
		return MatrixType::REAL_SYMMETRIC_INDEFINITE;
	}
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

bool StructuralMechanics3DSolverDataProvider::FETI::hasKernel(esint domain)
{
	if (_configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT) {
		return false;
	}
	return true;
}

int StructuralMechanics3DSolverDataProvider::FETI::initKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
//	for (esint d = 0; d < K.domains; ++d) {
//		if (K[d].type != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
//			eslog::error("Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC.\n");
//		}
//	}

	dnodes.resize(info::mesh->elements->ndomains);
	auto dmap = info::mesh->nodes->domains->cbegin();
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				dnodes[*d - info::mesh->elements->firstDomain].push_back(i);
			}
		}
	}

	N1.initDomains(K.domains);
	N2.initDomains(K.domains);
	RegMat.initDomains(K.domains);
	_RegMat = new MatrixCSRFETI();
	_RegMat->initDomains(info::mesh->elements->ndomains);

	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			N1[d].resize(K[d].nrows, 6);
		}
	}
	return 6;
}

void StructuralMechanics3DSolverDataProvider::FETI::fillKernels(MatrixCSRFETI &K, MatrixCSRFETI &M, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			std::vector<double> kk(K[d].nrows * K[d].nrows), mm(K[d].nrows * K[d].nrows);
			for (esint r = 0; r < K[d].nrows; ++r) {
				for (esint c = K[d].rows[r] - 1; c < K[d].rows[r + 1] - 1; ++c) {
					kk[r * K[d].ncols + K[d].cols[c] - 1] = K[d].vals[c];
					mm[r * M[d].ncols + M[d].cols[c] - 1] = M[d].vals[c];
				}
			}

			int m, info, itype = 1, il = 1, iu = 6, lwork = -1, iwork = -1;
			double abstol = 1e-6;
			std::vector<double> work(1), E(K[d].nrows), X(K[d].nrows * 6), res(K[d].nrows);
			std::vector<int> ifail(K[d].nrows);
			dsygvx(&itype, "V", "I", "L", &K[d].nrows, kk.data(), &K[d].nrows, mm.data(), &K[d].nrows, NULL, NULL, &il, &iu, &abstol, &m, E.data(), X.data(), &K[d].nrows, work.data(), &lwork, &iwork, ifail.data(), &info);
			work.resize(lwork = iwork = work.front());
			dsygvx(&itype, "V", "I", "L", &K[d].nrows, kk.data(), &K[d].nrows, mm.data(), &K[d].nrows, NULL, NULL, &il, &iu, &abstol, &m, E.data(), X.data(), &K[d].nrows, work.data(), &lwork, &iwork, ifail.data(), &info);

			N1[d].fillValues(X.data());
		}
	}
}

int StructuralMechanics3DSolverDataProvider::Hypre::numfnc()
{
	return 3;
}

void StructuralMechanics3DSolverDataProvider::Hypre::initKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N)
{
	N.initVectors(6);
	N.resize(K.nrows, K.nhalo, K.nneighbors);
}

void StructuralMechanics3DSolverDataProvider::Hypre::fillKernels(MatrixCSRDistributed &K, VectorsDenseDistributed &N)
{
	for (esint n = 0; n < info::mesh->nodes->size; n++) {
		Point p = info::mesh->nodes->coordinates->datatarray()[n];

		N[0][3 * n + 0] = 1;
		N[0][3 * n + 1] = 0;
		N[0][3 * n + 2] = 0;

		N[1][3 * n + 0] = 0;
		N[1][3 * n + 1] = 1;
		N[1][3 * n + 2] = 0;

		N[2][3 * n + 0] = 0;
		N[2][3 * n + 1] = 0;
		N[2][3 * n + 2] = 1;

		N[3][3 * n + 0] = -p.y;
		N[3][3 * n + 1] =  p.x;
		N[3][3 * n + 2] =    0;

		N[4][3 * n + 0] = -p.z;
		N[4][3 * n + 1] =    0;
		N[4][3 * n + 2] =  p.x;

		N[5][3 * n + 0] =    0;
		N[5][3 * n + 1] = -p.z;
		N[5][3 * n + 2] =  p.y;
	}
}
