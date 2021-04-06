
#include "structuralmechanics3d.provider.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/fetidatastore.h"
#include "math/matrix.type.h"
#include "math/matrix.dense.feti.h"
#include "math/matrix.csr.feti.h"
#include "math/matrix.csr.distributed.h"
#include "math/vector.dense.h"
#include "math/vector.dense.distributed.h"

#include <cmath>
#include <algorithm>

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

void StructuralMechanics3DSolverDataProvider::FETI::initKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	for (esint d = 0; d < K.domains; ++d) {
		if (K[d].type != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
			eslog::error("Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set FETI_REGULARIZATION = ALGEBRAIC.\n");
		}
	}

	size_t clusters = *std::max_element(info::mesh->elements->clusters.begin(), info::mesh->elements->clusters.end()) + 1;

	_cCenter = _cNorm = std::vector<Point>(clusters, Point(0, 0, 0));
	_cr44 = _cr45 = _cr46 = _cr55 = _cr56 = std::vector<double>(clusters, 0);
	_cNp = std::vector<size_t>(clusters, 0);

	_dCenter = _dNorm = std::vector<Point>(info::mesh->elements->ndomains, Point(0, 0, 0));
	_dr44 = _dr45 = _dr46 = _dr55 = _dr56 = std::vector<double>(info::mesh->elements->ndomains, 0);
	dnodes.resize(info::mesh->elements->ndomains);

	std::vector<double> cbuffer1(info::mesh->elements->ndomains, 0), cbuffer2(info::mesh->elements->ndomains, 0), cbuffer3(info::mesh->elements->ndomains, 0);

	auto dmap = info::mesh->nodes->domains->cbegin();
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				_dCenter[domain] += info::mesh->nodes->coordinates->datatarray()[i];
				dnodes[domain].push_back(i);
			}
		}
	}

	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		_cCenter[info::mesh->elements->clusters[d]] += _dCenter[d];
		_dCenter[d] = _dCenter[d] / dnodes[d].size();
		_cNp[info::mesh->elements->clusters[d]] += dnodes[d].size();
	}
	for (size_t c = 0; c < clusters; c++) {
		_cCenter[c] /= _cNp[c];
	}

	// Compute norm of column 4 (norm.x)
	dmap = info::mesh->nodes->domains->cbegin();
	std::vector<double> pnorm(info::mesh->elements->ndomains), pcnorm(info::mesh->elements->ndomains);
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				Point dp = info::mesh->nodes->coordinates->datatarray()[i] - _dCenter[domain];
				pnorm[domain] += dp.x * dp.x + dp.y * dp.y;
				Point cp = info::mesh->nodes->coordinates->datatarray()[i] - _cCenter[info::mesh->elements->clusters[domain]];
				pcnorm[domain] += cp.x * cp.x + cp.y * cp.y;
			}
		}
	}
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		_dNorm[d].x = std::sqrt(pnorm[d]);
		cbuffer1[d] += pcnorm[d];
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cNorm[info::mesh->elements->clusters[p]].x += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].x = std::sqrt(_cNorm[c].x);
	}

	// Compute coefficient r44, r45
	cbuffer1 = cbuffer2 = std::vector<double>(info::mesh->elements->ndomains, 0);
	dmap = info::mesh->nodes->domains->cbegin();
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				int cluster = info::mesh->elements->clusters[domain];

				Point dp = info::mesh->nodes->coordinates->datatarray()[i] - _dCenter[domain];
				_dr44[domain] += (-dp.y / _dNorm[domain].x) * (-dp.y / _dNorm[domain].x) + (dp.x / _dNorm[domain].x) * (dp.x / _dNorm[domain].x);
				_dr45[domain] += (-dp.y / _dNorm[domain].x) * (-dp.z);

				Point cp = info::mesh->nodes->coordinates->datatarray()[i] - _cCenter[cluster];
				cbuffer1[domain] += (-cp.y / _cNorm[cluster].x) * (-cp.y / _cNorm[cluster].x) + (cp.x / _cNorm[cluster].x) * (cp.x / _cNorm[cluster].x);
				cbuffer2[domain] += (-cp.y / _cNorm[cluster].x) * (-cp.z);
			}
		}
	}

	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cr44[info::mesh->elements->clusters[p]] += cbuffer1[p];
		_cr45[info::mesh->elements->clusters[p]] += cbuffer2[p];
	}

	// Compute norm of column 5 (norm.y)
	std::vector<double> dnorm(info::mesh->elements->ndomains), cnorm(info::mesh->elements->ndomains);
	cbuffer1 = std::vector<double>(info::mesh->elements->ndomains, 0);
	dmap = info::mesh->nodes->domains->cbegin();
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				int cluster = info::mesh->elements->clusters[domain];

				Point dp = info::mesh->nodes->coordinates->datatarray()[i] - _dCenter[domain];
				dnorm[domain] += (-dp.z - _dr45[domain] / _dr44[domain] * (-dp.y / _dNorm[domain].x)) * (-dp.z - _dr45[domain] / _dr44[domain] * (-dp.y / _dNorm[domain].x));
				dnorm[domain] += (    0 - _dr45[domain] / _dr44[domain] * ( dp.x / _dNorm[domain].x)) * (    0 - _dr45[domain] / _dr44[domain] * ( dp.x / _dNorm[domain].x));
				dnorm[domain] += dp.x * dp.x;

				Point cp = info::mesh->nodes->coordinates->datatarray()[i] - _cCenter[cluster];
				cnorm[domain] += (-cp.z - _cr45[cluster] / _cr44[cluster] * (-cp.y / _cNorm[cluster].x)) * (-cp.z - _cr45[cluster] / _cr44[cluster] * (-cp.y / _cNorm[cluster].x));
				cnorm[domain] += (    0 - _cr45[cluster] / _cr44[cluster] * ( cp.x / _cNorm[cluster].x)) * (    0 - _cr45[cluster] / _cr44[cluster] * ( cp.x / _cNorm[cluster].x));
				cnorm[domain] += cp.x * cp.x;
			}
		}
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_dNorm[p].y = std::sqrt(dnorm[p]);
		cbuffer1[p] = cnorm[p];
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cNorm[info::mesh->elements->clusters[p]].y += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].y = std::sqrt(_cNorm[c].y);
	}

	// Compute coefficient r46, r55, r56
	cbuffer1 = cbuffer2 = cbuffer3 = std::vector<double>(info::mesh->elements->ndomains, 0);
	std::vector<double> c5(info::mesh->elements->ndomains);
	dmap = info::mesh->nodes->domains->cbegin();
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				int cluster = info::mesh->elements->clusters[domain];
				double c5;

				Point dp = info::mesh->nodes->coordinates->datatarray()[i] - _dCenter[domain];
				_dr46[domain] += (dp.x / _dNorm[domain].x) * (-dp.z);
				c5 = (-dp.z - _dr45[domain] / _dr44[domain] * (-dp.y / _dNorm[domain].x)) / _dNorm[domain].y;
				_dr55[domain] += c5 * c5;
				_dr56[domain] += c5 * 0;
				c5 = (    0 - _dr45[domain] / _dr44[domain] * ( dp.x / _dNorm[domain].x)) / _dNorm[domain].y;
				_dr55[domain] += c5 * c5;
				_dr56[domain] += c5 * (-dp.z);
				c5 = ( dp.x -                                           0) / _dNorm[domain].y;
				_dr55[domain] += c5 * c5;
				_dr56[domain] += c5 * dp.y;

				Point cp = info::mesh->nodes->coordinates->datatarray()[i] - _cCenter[cluster];
				cbuffer1[domain] += (cp.x / _cNorm[cluster].x) * (-cp.z);
				c5 = (-cp.z - _cr45[cluster] / _cr44[cluster] * (-cp.y / _cNorm[cluster].x)) / _cNorm[cluster].y;
				cbuffer2[domain] += c5 * c5;
				cbuffer3[domain] += c5 * 0;
				c5 = (    0 - _cr45[cluster] / _cr44[cluster] * ( cp.x / _cNorm[cluster].x)) / _cNorm[cluster].y;
				cbuffer2[domain] += c5 * c5;
				cbuffer3[domain] += c5 * (-cp.z);
				c5 = ( cp.x -                                           0) / _cNorm[cluster].y;
				cbuffer2[domain] += c5 * c5;
				cbuffer3[domain] += c5 * cp.y;
			}
		}
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cr46[info::mesh->elements->clusters[p]] += cbuffer1[p];
		_cr55[info::mesh->elements->clusters[p]] += cbuffer2[p];
		_cr56[info::mesh->elements->clusters[p]] += cbuffer3[p];
	}

	// Compute norm of column 6 (norm.z)
	cbuffer1 = std::vector<double>(info::mesh->elements->ndomains, 0);
	std::vector<double> c6(info::mesh->elements->ndomains);
	dnorm.clear();
	cnorm.clear();
	dnorm.resize(info::mesh->elements->ndomains);
	cnorm.resize(info::mesh->elements->ndomains);
	dmap = info::mesh->nodes->domains->cbegin();
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++dmap) {
		for (auto d = dmap->begin(); d != dmap->end(); ++d) {
			if (info::mesh->elements->firstDomain <= *d && *d < info::mesh->elements->firstDomain + info::mesh->elements->ndomains) {
				esint domain = *d - info::mesh->elements->firstDomain;
				int cluster = info::mesh->elements->clusters[domain];
				double c6;

				Point dp = info::mesh->nodes->coordinates->datatarray()[i] - _dCenter[domain];
				c6 =     0 - _dr56[domain] / _dr55[domain] * (-dp.z - _dr45[domain] / _dr44[domain] * (-dp.y / _dNorm[domain].x)) / _dNorm[domain].y - _dr46[domain] / _dr44[domain] * (-dp.y / _dNorm[domain].x);
				dnorm[domain] += c6 * c6;
				c6 = -dp.z - _dr56[domain] / _dr55[domain] * (    0 - _dr45[domain] / _dr44[domain] * ( dp.x / _dNorm[domain].x)) / _dNorm[domain].y - _dr46[domain] / _dr44[domain] * ( dp.x / _dNorm[domain].x);
				dnorm[domain] += c6 * c6;
				c6 =  dp.y - _dr56[domain] / _dr55[domain] * ( dp.x -                                           0) / _dNorm[domain].y - _dr46[domain] / _dr44[domain] * (    0 / _dNorm[domain].x);
				dnorm[domain] += c6 * c6;

				Point cp = info::mesh->nodes->coordinates->datatarray()[i] - _cCenter[cluster];
				c6 =     0 - _cr56[cluster] / _cr55[cluster] * (-cp.z - _cr45[cluster] / _cr44[cluster] * (-cp.y / _cNorm[cluster].x)) / _cNorm[cluster].y - _cr46[cluster] / _cr44[cluster] * (-cp.y / _cNorm[cluster].x);
				cnorm[domain] += c6 * c6;
				c6 = -cp.z - _cr56[cluster] / _cr55[cluster] * (    0 - _cr45[cluster] / _cr44[cluster] * ( cp.x / _cNorm[cluster].x)) / _cNorm[cluster].y - _cr46[cluster] / _cr44[cluster] * ( cp.x / _cNorm[cluster].x);
				cnorm[domain] += c6 * c6;
				c6 =  cp.y - _cr56[cluster] / _cr55[cluster] * ( cp.x -                                           0) / _cNorm[cluster].y - _cr46[cluster] / _cr44[cluster] * (    0 / _cNorm[cluster].x);
				cnorm[domain] += c6 * c6;
			}
		}
	}

	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_dNorm[p].z = std::sqrt(dnorm[p]);
		cbuffer1[p] = cnorm[p];
	}
	for (esint p = 0; p < info::mesh->elements->ndomains; p++) {
		_cNorm[info::mesh->elements->clusters[p]].z += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].z = std::sqrt(_cNorm[c].z);
	}


	N1.initDomains(K.domains);
	N2.initDomains(K.domains);
	RegMat.initDomains(K.domains);
	_RegMat = new MatrixCSRFETI();
	_RegMat->initDomains(info::mesh->elements->ndomains);

	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			std::vector<esint> fixPoints;
		//	if (_BEMDomain[domain]) {
		//		fixPoints = std::vector<esint>(
		//				info::mesh->FETIData->surfaceFixPoints.begin() + info::mesh->FETIData->sFixPointsDistribution[domain],
		//				info::mesh->FETIData->surfaceFixPoints.begin() + info::mesh->FETIData->sFixPointsDistribution[domain + 1]);
		//	} else {
				fixPoints = std::vector<esint>(
						info::mesh->FETIData->innerFixPoints.begin() + info::mesh->FETIData->iFixPointsDistribution[d],
						info::mesh->FETIData->innerFixPoints.begin() + info::mesh->FETIData->iFixPointsDistribution[d + 1]);
		//	}

			MatrixCSR Nt(6, K[d].ncols, 9 * fixPoints.size());

			Nt.rows[0] = 1;
			Nt.rows[1] = Nt.rows[0] + fixPoints.size();
			Nt.rows[2] = Nt.rows[1] + fixPoints.size();
			Nt.rows[3] = Nt.rows[2] + fixPoints.size();
			Nt.rows[4] = Nt.rows[3] + 2 * fixPoints.size();
			Nt.rows[5] = Nt.rows[4] + 2 * fixPoints.size();
			Nt.rows[6] = Nt.rows[5] + 2 * fixPoints.size();

			auto n2DOF = [&] (esint node) {
				return std::lower_bound(dnodes[d].begin(), dnodes[d].end(), node) - dnodes[d].begin();
			};

			esint cindex = 0;
			for (size_t c = 0; c < 3; c++) {
				for (size_t i = 0; i < fixPoints.size(); i++, cindex++) {
					Nt.cols[cindex] = 3 * n2DOF(fixPoints[i]) + c + 1;
					Nt.vals[cindex] = 1;
				}
			}

			for (size_t i = 0; i < fixPoints.size(); i++, cindex += 2) {
				const Point &p = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
				Nt.cols[cindex]     = 3 * n2DOF(fixPoints[i]) + 0 + 1;
				Nt.cols[cindex + 1] = 3 * n2DOF(fixPoints[i]) + 1 + 1;
				Nt.vals[cindex]     = -p.y;
				Nt.vals[cindex + 1] =  p.x;
			}

			for (size_t i = 0; i < fixPoints.size(); i++, cindex += 2) {
				const Point &p = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
				Nt.cols[cindex]     = 3 * n2DOF(fixPoints[i]) + 0 + 1;
				Nt.cols[cindex + 1] = 3 * n2DOF(fixPoints[i]) + 2 + 1;
				Nt.vals[cindex]     = -p.z;
				Nt.vals[cindex + 1] =  p.x;
			}

			for (size_t i = 0; i < fixPoints.size(); i++, cindex += 2) {
				const Point &p = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]];
				Nt.cols[cindex]     = 3 * n2DOF(fixPoints[i]) + 1 + 1;
				Nt.cols[cindex + 1] = 3 * n2DOF(fixPoints[i]) + 2 + 1;
				Nt.vals[cindex]     = -p.z;
				Nt.vals[cindex + 1] =  p.y;
			}

			// N * (tran(N) * N)^-1 * tran(N)
			//
			// AX = B => X = A^-1B
			//
			// if A = tran(N) * N;
			// then X = (tran(N) * N)^-1 * N
			// RegMat = tran(N) * X

			MatrixCSR N;
			Nt.transposeTo(&N);
			MatrixCSR A;

			A.multiply(Nt, N);
			A.removeLower(MatrixType::REAL_SYMMETRIC_INDEFINITE);
			MatrixDense _X(N.nrows, N.ncols), B = N;
			A.solve(B, _X);
			MatrixCSR X = _X;
			X.transpose();
			_RegMat->at(d)->multiply(N, X);
			_RegMat->at(d)->removeLower(MatrixType::REAL_SYMMETRIC_INDEFINITE);

			RegMat[d].shallowCopyStructure(_RegMat->at(d));

			N1[d].resize(K[d].nrows, 6);
		}
	}
}

void StructuralMechanics3DSolverDataProvider::FETI::fillKernels(MatrixCSRFETI &K, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			VectorDense diag(K[d].nrows, N1[d].vals);
			K[d].fillDiagonal(&diag);

			RegMat[d].fillData(_RegMat->at(d));
			RegMat[d].scale(diag.max());

			Point center = _dCenter[d], norm = _dNorm[d];
			double r44 = _dr44[d], r45 = _dr45[d], r46 = _dr46[d], r55 = _dr55[d], r56 = _dr56[d];
			size_t np = dnodes[d].size();

			if (ortogonalizeCluster) {
				size_t cluster = info::mesh->elements->clusters[d];
				center = _cCenter[cluster], norm = _cNorm[cluster];
				r44 = _cr44[cluster], r45 = _cr45[cluster], r46 = _cr46[cluster], r55 = _cr55[cluster], r56 = _cr56[cluster];
				np = _cNp[cluster];
			} else {
				center = _dCenter[d], norm = _dNorm[d];
				r44 = _dr44[d], r45 = _dr45[d], r46 = _dr46[d], r55 = _dr55[d], r56 = _dr56[d];
				np = dnodes[d].size();
			}

			double v = 1 / std::sqrt(np);
			for (size_t i = 0; i < dnodes[d].size(); i++) {
				Point p = info::mesh->nodes->coordinates->datatarray()[dnodes[d][i]] - center;

				N1[d][3 * i + 0][0] = v;
				N1[d][3 * i + 1][0] = 0;
				N1[d][3 * i + 2][0] = 0;

				N1[d][3 * i + 0][1] = 0;
				N1[d][3 * i + 1][1] = v;
				N1[d][3 * i + 2][1] = 0;

				N1[d][3 * i + 0][2] = 0;
				N1[d][3 * i + 1][2] = 0;
				N1[d][3 * i + 2][2] = v;

				N1[d][3 * i + 0][3] = -p.y / norm.x;
				N1[d][3 * i + 1][3] =  p.x / norm.x;
				N1[d][3 * i + 2][3] =             0;

				N1[d][3 * i + 0][4] = (-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y;
				N1[d][3 * i + 1][4] = (   0 - r45 / r44 * ( p.x / norm.x)) / norm.y;
				N1[d][3 * i + 2][4] = ( p.x - r45 / r44 * (   0 / norm.x)) / norm.y;

				N1[d][3 * i + 0][5] = (   0 - r56 / r55 * ((-p.z - r45 / r44 * (-p.y / norm.x)) / norm.y) - r46 / r44 * (-p.y / norm.x)) / norm.z;
				N1[d][3 * i + 1][5] = (-p.z - r56 / r55 * ((   0 - r45 / r44 * ( p.x / norm.x)) / norm.y) - r46 / r44 * ( p.x / norm.x)) / norm.z;
				N1[d][3 * i + 2][5] = ( p.y - r56 / r55 * (( p.x - r45 / r44 * (   0 / norm.x)) / norm.y) - r46 / r44 * (   0 / norm.x)) / norm.z;
			}
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
