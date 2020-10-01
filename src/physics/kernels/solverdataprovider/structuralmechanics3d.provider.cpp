
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

	size_t num_directions = _configuration.feti.num_directions;

	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			N1[d].resize(2 * K[d].nrows, 6 * num_directions);
		}
	}
	return 6;
}

void StructuralMechanics3DSolverDataProvider::FETI::fillKernels(MatrixCSRFETI &K, MatrixCSRFETI &M, MatrixDenseFETI &N1, MatrixDenseFETI &N2, MatrixCSRFETI &RegMat, bool ortogonalizeCluster)
{
	printf("fill kernels\n");

	size_t num_directions = _configuration.feti.num_directions;
/*
	std::vector<Point> wave_directions {
        Point{1, 0, 0},
        Point{0, 1, 0},
        Point{0, 0, 1},
        Point{7.071067811865475e-01, 7.071067811865475e-01, 0.000000000000000e+00},
        Point{7.071067811865475e-01, 0.000000000000000e+00, 7.071067811865475e-01},
        Point{0.000000000000000e+00, 7.071067811865475e-01, 7.071067811865475e-01},
        Point{5.773502691896258e-01, 5.773502691896258e-01, 5.773502691896258e-01}
	};

	std::vector<Point> As1 {
        Point{0, 1, 0},
        Point{-1, 0, 0},
        Point{0, 1, 0},
        Point{-7.071067811865475e-01, 7.071067811865476e-01, 0.000000000000000e+00},
        Point{0, 1, 0},
        Point{-7.071067811865476e-01, 4.999999999999999e-01, -5.000000000000001e-01},
        Point{-5.773502691896258e-01, 7.886751345948129e-01, -2.113248654051872e-01}
	};

	std::vector<Point> As2 {
        Point{0, 0, 1},
        Point{0, 0, 1},
        Point{-1, 0, 0},
        Point{0, 0, 1},
        Point{-7.071067811865475e-01, 0.000000000000000e+00, 7.071067811865476e-01},
        Point{-7.071067811865476e-01, -5.000000000000001e-01, 4.999999999999999e-01},
        Point{-5.773502691896258e-01, -2.113248654051872e-01, 7.886751345948129e-01}
	};
*/
    std::vector<Point> wave_directions {
        Point{1.000000, 0.000000, 0.000000},
        Point{0.000000, 1.000000, 0.000000},
        Point{0.000000, 0.000000, 1.000000},
        Point{1.000000, 1.000000, 0.000000},
        Point{1.000000, -1.000000, 0.000000},
        Point{1.000000, 1.000000, 1.000000},
        Point{-1.000000, 1.000000, 1.000000},
        Point{1.000000, -1.000000, 1.000000},
        Point{-1.000000, -1.000000, 1.000000},
        Point{1.000000, 0.000000, 1.000000},
        Point{-1.000000, 0.000000, 1.000000},
        Point{0.000000, 1.000000, 1.000000},
        Point{0.000000, -1.000000, 1.000000}
    };

    std::vector<Point> As1 {
        Point{0.000000, 1.000000, 0.000000},
        Point{-1.000000, 0.000000, 0.000000},
        Point{0.000000, 1.000000, 0.000000},
        Point{-0.707107, 0.707107, 0.000000},
        Point{0.707107, 0.707107, 0.000000},
        Point{-0.577350, 0.788675, -0.211325},
        Point{0.577350, 0.788675, -0.211325},
        Point{0.577350, 0.788675, 0.211325},
        Point{-0.577350, 0.788675, 0.211325},
        Point{0.000000, 1.000000, 0.000000},
        Point{0.000000, 1.000000, 0.000000},
        Point{-0.707107, 0.500000, -0.500000},
        Point{0.707107, 0.500000, 0.500000}
    };

    std::vector<Point> As2 {
        Point{0.000000, 0.000000, 1.000000},
        Point{0.000000, 0.000000, 1.000000},
        Point{-1.000000, 0.000000, 0.000000},
        Point{0.000000, 0.000000, 1.000000},
        Point{0.000000, 0.000000, 1.000000},
        Point{-0.577350, -0.211325, 0.788675},
        Point{0.577350, -0.211325, 0.788675},
        Point{-0.577350, 0.211325, 0.788675},
        Point{0.577350, 0.211325, 0.788675},
        Point{-0.707107, 0.000000, 0.707107},
        Point{0.707107, 0.000000, 0.707107},
        Point{-0.707107, -0.500000, 0.500000},
        Point{-0.707107, 0.500000, 0.500000}
    };

    double omega = 2*3.141592*100;
//    std::FILE *fp = std::fopen("omega.txt", "r");
//    std::fscanf(fp, "%lf", &omega);
//    std::fclose(fp);
    std::printf("omega = %f\n", omega);
    //const double omega = std::sqrt(100.0) / (2*3.1415926);
    double rho = 7850.0;
    double nu = 0.3;
    double E = 2e11;
    // const double Lambda = (mu * (E - 2 * mu)) / (3 * mu - E);
    // const double Lambda = (E * mu) / ((1.0 + mu)*(1.0 - 2.0 * mu));
    //mu = E / (2.0*(1+mu));
    double Lambda = E*nu / ((1.0+nu)*(1.0-2.0*nu));
//    double mu = E*(1.0-nu) / ((1.0+nu)*(1.0-2.0*nu));
    double mu = E / (2*(1.0+nu));

    double kp = std::sqrt((rho * omega * omega) / (Lambda + 2.0 * mu));
    double ks = std::sqrt((rho * omega * omega) / mu);

    /*
	std::printf("omega = %f\n", omega);
	std::printf("rho = %f\n", rho);
	std::printf("E = %f\n", E);
    */
    std::printf("mu = %f\n", mu);
	std::printf("Lambda = %f\n", Lambda);

	#pragma omp parallel for
	for (esint d = 0; d < K.domains; ++d) {
		if (hasKernel(d)) {
			for (size_t i = 0; i < dnodes[d].size(); i++) {
				Point p = info::mesh->nodes->coordinates->datatarray()[dnodes[d][i]];

				for (size_t j = 0; j < num_directions; ++j) {
					double tp = kp * (p * wave_directions[j]);
					double ts = ks * (p * wave_directions[j]);

					Point a0  = wave_directions[j];
					Point s1 = As1[j];
					Point s2 = As2[j];

					// Re x
					N1[d][6 * i + 0][6 * j + 0] = a0.x * std::cos(tp);
					N1[d][6 * i + 0][6 * j + 1] = s1.x * std::cos(ts);
					N1[d][6 * i + 0][6 * j + 2] = s2.x * std::cos(ts);

					// - Im x
					N1[d][6 * i + 0][6 * j + 3] = -a0.x * std::sin(tp);
					N1[d][6 * i + 0][6 * j + 4] = -s1.x * std::sin(ts);
					N1[d][6 * i + 0][6 * j + 5] = -s2.x * std::sin(ts);


					// Re y
					N1[d][6 * i + 1][6 * j + 0] = a0.y * std::cos(tp);
					N1[d][6 * i + 1][6 * j + 1] = s1.y * std::cos(ts);
					N1[d][6 * i + 1][6 * j + 2] = s2.y * std::cos(ts);

					// - Im y
					N1[d][6 * i + 1][6 * j + 3] = -a0.y * std::sin(tp);
					N1[d][6 * i + 1][6 * j + 4] = -s1.y * std::sin(ts);
					N1[d][6 * i + 1][6 * j + 5] = -s2.y * std::sin(ts);


					// Re z
					N1[d][6 * i + 2][6 * j + 0] = a0.z * std::cos(tp);
					N1[d][6 * i + 2][6 * j + 1] = s1.z * std::cos(ts);
					N1[d][6 * i + 2][6 * j + 2] = s2.z * std::cos(ts);

					// - Im z
					N1[d][6 * i + 2][6 * j + 3] = -a0.z * std::sin(tp);
					N1[d][6 * i + 2][6 * j + 4] = -s1.z * std::sin(ts);
					N1[d][6 * i + 2][6 * j + 5] = -s2.z * std::sin(ts);


					// Im x
					N1[d][6 * i + 3][6 * j + 0] = a0.x * std::sin(tp);
					N1[d][6 * i + 3][6 * j + 1] = s1.x * std::sin(ts);
					N1[d][6 * i + 3][6 * j + 2] = s2.x * std::sin(ts);

					// Re x
					N1[d][6 * i + 3][6 * j + 3] = a0.x * std::cos(tp);
					N1[d][6 * i + 3][6 * j + 4] = s1.x * std::cos(ts);
					N1[d][6 * i + 3][6 * j + 5] = s2.x * std::cos(ts);


					// Im y
					N1[d][6 * i + 4][6 * j + 0] = a0.y * std::sin(tp);
					N1[d][6 * i + 4][6 * j + 1] = s1.y * std::sin(ts);
					N1[d][6 * i + 4][6 * j + 2] = s2.y * std::sin(ts);

					// Re y
					N1[d][6 * i + 4][6 * j + 3] = a0.y * std::cos(tp);
					N1[d][6 * i + 4][6 * j + 4] = s1.y * std::cos(ts);
					N1[d][6 * i + 4][6 * j + 5] = s2.y * std::cos(ts);


					// Im z
					N1[d][6 * i + 5][6 * j + 0] = a0.z * std::sin(tp);
					N1[d][6 * i + 5][6 * j + 1] = s1.z * std::sin(ts);
					N1[d][6 * i + 5][6 * j + 2] = s2.z * std::sin(ts);

					// Re z
					N1[d][6 * i + 5][6 * j + 3] = a0.z * std::cos(tp);
					N1[d][6 * i + 5][6 * j + 4] = s1.z * std::cos(ts);
					N1[d][6 * i + 5][6 * j + 5] = s2.z * std::cos(ts);
				}
			}
		}
	}
}

std::vector<Point> StructuralMechanics3DSolverDataProvider::FETI::getWaveDirections(size_t dir_steps)
{
	// number of points on the surface of the cube with a given number of steps
	size_t num_direction_points = 6 * dir_steps * (dir_steps - 2) + 8;

	std::vector<Point> wave_directions;
	wave_directions.reserve(num_direction_points);

	const double step_size = 2.0 / (dir_steps - 1);

	for (size_t i = 0; i < dir_steps; ++i) {
		for (size_t j = 0; j < dir_steps; ++j) {
			// top
			wave_directions.push_back(Point(-1.0 + step_size * i, -1.0 + step_size * j, 1.0));

			// bottom
			wave_directions.push_back(Point(-1.0 + step_size * i, -1.0 + step_size * j, -1.0));
		}
	}

	// sides
	for (size_t i = 1; i < dir_steps-1; ++i) {
		for (size_t j = 0; j < dir_steps-1; ++j) {
			wave_directions.push_back(Point(-1.0 + step_size * j, -1.0, -1.0 + step_size * i));
			wave_directions.push_back(Point(1.0, -1.0 + step_size * j, -1.0 + step_size * i));
			wave_directions.push_back(Point(1.0 - step_size * j, 1.0, -1.0 + step_size * i));
			wave_directions.push_back(Point(-1.0, 1.0 - step_size * j, -1.0 + step_size * i));
		}
	}

	// normalize all vectors to unit length
	for (size_t i = 0; i < num_direction_points; ++i) {
		wave_directions[i].normalize();
	}

	return wave_directions;
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
