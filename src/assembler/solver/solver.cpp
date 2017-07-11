
#include "solver.h"

#include "../../basis/logging/constants.h"
#include "../../basis/utilities/utils.h"
#include "../../configuration/environment.h"
#include "../../configuration/output.h"
#include "../../configuration/solver/espresooptions.h"
#include "../physics/physics.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"

#include "../../mesh/structures/mesh.h"
#include "../../output/resultstore.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../solver/generic/LinearSolver.h"

using namespace espreso;

Solver::Solver(
		const std::string &name,
		Mesh *mesh,
		Physics* physics,
		LinearSolver* linearSolver,
		output::Store* store,
		double duration,
		Matrices restriction)
: SolverBase(name, physics->name(), mesh, duration), physics(physics), instance(physics->instance()), linearSolver(linearSolver), _store(store), _restriction(~restriction)
{

}

Solver::~Solver()
{

}

static std::string mNames(espreso::Matrices matrices)
{
	return
	std::string(matrices & espreso::Matrices::K      ? "K "      : "") +
	std::string(matrices & espreso::Matrices::M      ? "M "      : "") +
	std::string(matrices & espreso::Matrices::R      ? "R "      : "") +
	std::string(matrices & espreso::Matrices::f      ? "f "      : "") +
	std::string(matrices & espreso::Matrices::B0     ? "B0 "     : "") +
	std::string(matrices & espreso::Matrices::B1     ? "B1 "     : "") +
	std::string(matrices & espreso::Matrices::B1c    ? "B1c "    : "") +
	std::string(matrices & espreso::Matrices::primal ? "Primar " : "") +
	std::string(matrices & espreso::Matrices::dual   ? "Dual "   : "");
}

void Solver::storeData(const Step &step, std::vector<SparseMatrix> &matrices, const std::string &name, const std::string &description)
{
	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing " << description;
		for (size_t d = 0; d < matrices.size(); d++) {
			std::ofstream os(Logging::prepareFile(d, name));
			os << matrices[d];
			os.close();
		}
	}
}

void Solver::storeData(const Step &step, std::vector<std::vector<double> > &vectors, const std::string &name, const std::string &description)
{
	if (environment->print_matrices) {
		ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing " << description;
		for (size_t d = 0; d < vectors.size(); d++) {
			std::ofstream os(Logging::prepareFile(d, name));
			os << vectors[d];
			os.close();
		}
	}
}

void Solver::storeSolution(const Step &step)
{
	std::vector<Solution*> solutions;
	for (size_t i = 0; i < physics->solutions().size(); i++) {
		solutions.push_back(instance->solutions[physics->solutions()[i]]);
	}

	_store->storeSolution(step, solutions, physics->properties());
}

void Solver::storeSubSolution(const Step &step)
{
	if (_store->configuration().iterations) {
		_store->storeSolution(step, instance->solutions);
	}
}

void Solver::preprocessData(const Step &step)
{
	physics->preprocessData(step);
}

void Solver::updateMatrices(const Step &step, Matrices matrices)
{
	updateMatrices(step, matrices, instance->solutions);
}

void Solver::updateMatrices(const Step &step, Matrices matrices, const std::vector<Solution*> &solution)
{
	matrices &= _restriction;
	if (matrices & ~(Matrices::K | Matrices::M | Matrices::R | Matrices::f)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid matrices passed to Solver::assembleMatrices.";
	}
	if (!matrices) {
		return;
	}

	ESINFO(PROGRESS2) << physics->name() << (solution.size() ? " updates" : " assembles") << " matrices: " << mNames(matrices);

	TimeEvent time(std::string(solution.size() ? "Updates" : "Assembles") + " matrices " + mNames(matrices) + " by " + physics->name()); time.start();
	physics->updateMatrix(step, matrices, solution);
	storeData(step, physics->instance()->K, "origK", "original stiffness matrices K");
	time.endWithBarrier(); _timeStatistics->addEvent(time);
}

void Solver::subtractDirichlet()
{
	TimeEvent time(std::string("Subtract dirichlet")); time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < physics->instance()->domains; d++) {
		for (size_t j = 0; j < physics->instance()->B1[d].J_col_indices.size(); j++) {
			if (physics->instance()->B1[d].I_row_indices[j] > (eslocal)physics->instance()->block[Instance::CONSTRAINT::DIRICHLET]) {
				break;
			}
			physics->instance()->B1c[d][j] -= physics->instance()->primalSolution[d][physics->instance()->B1[d].J_col_indices[j] - 1];
		}
	}
	time.endWithBarrier(); _timeStatistics->addEvent(time);
}

void Solver::sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::string &description)
{
	std::vector<size_t> prefix(x.size());
	for (size_t i = 0; i < x.size(); i++) {
		prefix[i] = x[i].size();
	}
	sum(z, a, x, b, y, prefix, description);
}

void Solver::sum(std::vector<std::vector<double> > &z, double a, const std::vector<std::vector<double> > &x, double b, const std::vector<std::vector<double> > &y, const std::vector<size_t> &prefix, const std::string &description)
{
	ESINFO(PROGRESS2) << "Compute " << description;
	TimeEvent time("Compute " + description); time.start();
	if (z.size() == 0) {
		z.resize(x.size());
	}
	#pragma omp parallel for
	for (size_t d = 0; d < x.size(); d++) {
		if (z[d].size() == 0) {
			z[d].resize(x[d].size());
		}
		if (x[d].size() != y[d].size() || z[d].size() != x[d].size()) {
			ESINFO(ERROR) << "ESPRESO internal error while " << description << ". Vectors have different dimension.";
		}
		for (size_t i = 0; i < x[d].size() && i < prefix[d]; i++) {
			z[d][i] = a * x[d][i] + b * y[d][i];
		}
	}
	time.endWithBarrier(); _timeStatistics->addEvent(time);
}

void Solver::sum(std::vector<SparseMatrix> &A, double beta, std::vector<SparseMatrix> &B, const std::string &description)
{
	ESINFO(PROGRESS2) << "Compute " << description;
	TimeEvent time("Compute " + description); time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < physics->instance()->domains; d++) {
		A[d].MatAddInPlace(B[d], 'N', beta);
	}
	time.endWithBarrier(); _timeStatistics->addEvent(time);
}

void Solver::multiply(std::vector<std::vector<double> > &y, std::vector<SparseMatrix> &A, std::vector<std::vector<double> > &x, const std::string &description)
{
	ESINFO(PROGRESS2) << "Compute " << description;
	TimeEvent time("Compute " + description); time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < physics->instance()->domains; d++) {
		A[d].MatVec(x[d], y[d], 'N', 0, 0, 0);
	}
	time.endWithBarrier(); _timeStatistics->addEvent(time);
}

void Solver::multiply(std::vector<std::vector<double> > &x, double a, const std::string &description)
{
	ESINFO(PROGRESS2) << "Compute " << description;
	TimeEvent time("Compute " + description); time.start();
	#pragma omp parallel for
	for (size_t d = 0; d < x.size(); d++) {
		for (size_t i = 0; i < x[d].size(); i++) {
			x[d][i] *= a;
		}
	}
	time.endWithBarrier(); _timeStatistics->addEvent(time);
}

void Solver::regularizeMatrices(const Step &step, Matrices matrices)
{
	if (matrices & ~(Matrices::K)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid matrices passed to Solver::regularizeMatrices.";
	}

	instance->computeKernelsCallback = [&] (REGULARIZATION regularization, size_t scSize) {

		ESINFO(PROGRESS2) << physics->name() << " regularizes matrices: " << mNames(Matrices::K);

		TimeEvent time("Regularization of matrices " + mNames(Matrices::K) + " by " + physics->name()); time.start();
		physics->makeStiffnessMatricesRegular(regularization, scSize);
		time.endWithBarrier(); _timeStatistics->addEvent(time);

		storeData(step, physics->instance()->N1, "N1", "N1");
		storeData(step, physics->instance()->N2, "N2", "N2");
		storeData(step, physics->instance()->RegMat, "RegMat", "RegMat");
	};

	instance->computeKernelCallback = [&] (REGULARIZATION regularization, size_t scSize, size_t domain) {
		physics->makeStiffnessMatrixRegular(regularization, scSize, domain);
	};
}

void Solver::composeGluing(const Step &step, Matrices matrices)
{
	if (matrices & ~(Matrices::B0 | Matrices::B1)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid matrices passed to Solver::composeGluing.";
	}
	if (!matrices) {
		return;
	}

	instance->assembleB0Callback = [&] (B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
		ESINFO(PROGRESS2) << "Compose gluing matrices " << mNames( Matrices::B0) << "for " << physics->name();
		TimeEvent time("Composition of gluing matrices" + mNames(Matrices::B0) + "for " + physics->name()); time.start();
		physics->instance()->B0.clear();
		physics->instance()->B0.resize(physics->instance()->domains);
		for (size_t d = 0; d < physics->instance()->domains; d++) {
			physics->instance()->B0[d].type = 'G';
			physics->instance()->B0subdomainsMap[d].clear();
		}
		switch (type) {
		case B0_TYPE::CORNERS:
			physics->assembleB0FromCorners(step);
			break;
		case B0_TYPE::KERNELS:
			physics->assembleB0FromKernels(step, kernels);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
		}
		time.endWithBarrier(); _timeStatistics->addEvent(time);

		storeData(step, physics->instance()->B0, "B0", "B0");
	};

	if (matrices & Matrices::B1) {
		ESINFO(PROGRESS2) << "Compose gluing matrices " << mNames( Matrices::B1) << "for " << physics->name();
		TimeEvent time("Composition of gluing matrices" + mNames(Matrices::B1) + "for " + physics->name()); time.start();
		// TODO: create update method
		physics->instance()->B1.clear();
		physics->instance()->B1.resize(physics->instance()->domains);
		physics->instance()->inequality.clear();
		physics->instance()->inequality.resize(physics->instance()->domains);
		physics->instance()->B1clustersMap.clear();
		for (size_t d = 0; d < physics->instance()->domains; d++) {
			physics->instance()->B1[d].type = 'G';
			physics->instance()->B1c[d].clear();
			physics->instance()->LB[d].clear();
			physics->instance()->B1duplicity[d].clear();
			physics->instance()->inequalityC[d].clear();
			physics->instance()->B1subdomainsMap[d].clear();
		}
		physics->instance()->block.clear();
		physics->instance()->block.resize(3, 0);
		physics->assembleB1(step, linearSolver->configuration.redundant_lagrange, linearSolver->configuration.scaling);
		time.endWithBarrier(); _timeStatistics->addEvent(time);
	}

}

void Solver::processSolution(const Step &step)
{
	ESINFO(PROGRESS2) << "Run " << physics->name() << " post-processing";
	TimeEvent time(physics->name() + " post-processing"); time.start();
	physics->processSolution(step);
	time.endWithBarrier(); _timeStatistics->addEvent(time);

	storeData(step, physics->instance()->primalSolution, "solution", "solution");
}

void Solver::initLinearSolver(const Step &step)
{
	storeData(step, physics->instance()->K, "K", "stiffness matrices K");
	storeData(step, physics->instance()->M, "M", "mass matrices M");
	storeData(step, physics->instance()->R, "R", "residual forces R");
	storeData(step, physics->instance()->f, "f", "right-hand side");

	storeData(step, physics->instance()->B1, "B1", "B1");
	storeData(step, physics->instance()->B1c, "B1c", "B1c");
	storeData(step, physics->instance()->B1duplicity, "B1duplicity", "B1duplicity");

	ESINFO(PROGRESS2) << "Initialization of linear solver";

	TimeEvent timeSolver("Initialize linear solver"); timeSolver.startWithBarrier();
	linearSolver->init();
	timeSolver.end(); _timeStatistics->addEvent(timeSolver);
}

void Solver::updateLinearSolver(const Step &step, Matrices matrices)
{
	storeData(step, physics->instance()->K, "K", "stiffness matrices K");
	storeData(step, physics->instance()->M, "M", "mass matrices M");
	storeData(step, physics->instance()->R, "R", "residual forces R");
	storeData(step, physics->instance()->f, "f", "right-hand side");

	storeData(step, physics->instance()->B1, "B1", "B1");
	storeData(step, physics->instance()->B1c, "B1c", "B1c");
	storeData(step, physics->instance()->B1duplicity, "B1duplicity", "B1duplicity");

	ESINFO(PROGRESS2) << "Updating of linear solver";

	TimeEvent timeSolver("Update linear solver"); timeSolver.startWithBarrier();
	linearSolver->update(matrices);
	timeSolver.end(); _timeStatistics->addEvent(timeSolver);
}

void Solver::runLinearSolver(const Step &step)
{
	if (_store->configuration().FETI_data) {
		_store->storeFETIData(step, *instance);
	}

	ESINFO(PROGRESS2) << "Solve system";

	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	linearSolver->run();
	timeSolve.endWithBarrier(); _timeStatistics->addEvent(timeSolve);
}

void Solver::finalizeLinearSolver(const Step &step)
{
	ESINFO(PROGRESS2) << "Finalize " << _name << " solver";

	linearSolver->finilize();

	_store->finalize();

	_timeStatistics->totalTime.endWithBarrier();
	_timeStatistics->printStatsMPI();
}

double Solver::lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext, Physics *physics, const Step &step)
{
	auto multiply = [] (const std::vector<std::vector<double> > &v1, const std::vector<std::vector<double> > &v2) {
		double cmul = 0, gmul;

		#pragma omp parallel for
		for (size_t d = 0; d < v1.size(); d++) {
			double dmul = 0;
			for (size_t i = 0; i < v1[d].size(); i++) {
				dmul += v1[d][i] * v2[d][i];
			}
			#pragma omp atomic
			cmul += dmul;
		}
		MPI_Allreduce(&cmul, &gmul, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
		return gmul;
	};

	double a = 0, b = 1, alpha = 1;
	double fa = 0, fb = 0, fx = 0, faStart = 0;

	std::vector<std::vector<double> > solution = deltaU;
	std::vector<std::vector<double> > F_ext_r = F_ext;

	for (size_t i = 0; i < 6; i++) {
		sum(solution, 1, U, alpha, deltaU, "u = u + " + ASCII::alpha + " " + ASCII::DELTA + "u (line search)");

		solution.swap(physics->instance()->primalSolution);
		physics->updateMatrix(step, Matrices::R, physics->instance()->solutions);
		solution.swap(physics->instance()->primalSolution);

		if (i == 0) {
			faStart = multiply(deltaU, physics->instance()->f);
			sum(F_ext_r, 1, F_ext, -1, instance->R, "F_ext - R");

			fb = multiply(deltaU, F_ext_r);
			if ((faStart < 0 && fb < 0) || (faStart >= 0 && fb >= 0)) {
				return alpha;
			}
			fa = faStart;
		} else {
			sum(F_ext_r, 1, F_ext, -1, instance->R, "F_ext - R");
			fx = multiply(deltaU, F_ext_r);
			if (fa * fx < 0) {
				b = alpha;
				fb = fx;
			} else if (fb * fx < 0) {
				a = alpha;
				fa = fx;
			}

			if (fabs(fx) <= 0.5 * faStart) {
				alpha = a - fa * ((b - a ) / (fb - fa));
				break;
			}
		}

		alpha = a - fa * ((b - a ) / (fb - fa));
	}

	if (alpha < 0.1) {
		alpha = 0.1;
	}
	if (alpha > .99) {
		alpha = 1;
	}

	sum(solution, 0, U, alpha, deltaU, ASCII::DELTA + "u = " + ASCII::alpha + " * " + ASCII::DELTA + "u (line search)");
	solution.swap(deltaU);
	return alpha;
}

double Solver::maxAbsValue(const std::vector<std::vector<double> > &v) const
{
	double max = 0;
	for (size_t p = 0; p < v.size(); p++) {
		max = std::max(max, std::fabs(*std::max_element(v[p].begin(), v[p].end(), [] (const double v1, const double v2) { return std::fabs(v1) < std::fabs(v2); })));
	}

	double gmax;
	MPI_Allreduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);
	return gmax;
}

