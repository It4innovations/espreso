
#include "solver.h"

#include "../../basis/utilities/utils.h"
#include "../../configuration/environment.h"
#include "../../configuration/output.h"
#include "../../configuration/solver/espresooptions.h"
#include "../physics/physics.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"

#include "../../mesh/structures/mesh.h"

#include "../../output/vtk/vtk.h"

#include "../../solver/generic/SparseMatrix.h"
#include "../../solver/generic/LinearSolver.h"

using namespace espreso;

Solver::Solver(
		const std::string &name,
		Mesh *mesh,
		Physics* physics,
		LinearSolver* linearSolver,
		store::ResultStore* store,
		Matrices restriction)
: physics(physics), linearSolver(linearSolver), _name(name), _mesh(mesh), _store(store), _restriction(~restriction)
{
	_timeStatistics = new TimeEval(physics->name() + " solved by " + _name + " solver overall timing");
	_timeStatistics->totalTime.startWithBarrier();
}

Solver::~Solver()
{
	delete _timeStatistics;
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
	std::string(matrices & espreso::Matrices::primar ? "Primar " : "") +
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
	std::stringstream ss;
	ss << "_s" << step.step << "i" << step.iteration;
	for (size_t s = 0; s < physics->instance()->solutions.size(); s++) {
		_store->storeValues(physics->instance()->solutions[s]->name + ss.str(), physics->instance()->solutions[s]->properties, physics->instance()->solutions[s]->data, physics->instance()->solutions[s]->eType);
	}
}

void Solver::storeSubSolution(const Step &step)
{
	if (_store->configuration().substeps) {
		std::stringstream ss;
		ss << "_s" << step.step << "i" << step.iteration << "_" << step.substep;
		for (size_t s = 0; s < physics->instance()->solutions.size(); s++) {
			_store->storeValues(physics->instance()->solutions[s]->name + ss.str(), physics->instance()->solutions[s]->properties, physics->instance()->solutions[s]->data, physics->instance()->solutions[s]->eType);
		}
	}
}

void Solver::assembleMatrices(const Step &step, Matrices matrices)
{
	updateMatrices(step, matrices, {});
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
	time.endWithBarrier(); _timeStatistics->addEvent(time);

	if (matrices & Matrices::K) {
		storeData(step, physics->instance()->K, "K", "stiffness matrices K");
	}
	if (matrices & Matrices::M) {
		storeData(step, physics->instance()->M, "M", "mass matrices K");
	}
	if (matrices & Matrices::R) {
		storeData(step, physics->instance()->R, "R", "residual forces R");
	}
	if (matrices & Matrices::f) {
		storeData(step, physics->instance()->f, "f", "right-hand side");
	}
}

void Solver::updateVector(const Step &step, Matrices v1, Matrices v2, double alpha, double beta)
{
	ESINFO(PROGRESS2) << "Update vector: " << mNames(v1) << "= " << alpha << " * " << mNames(v1) << (beta < 0 ? "- " : "+ ") << fabs(beta) << " * " << mNames(v2);
	if (v1 & Matrices::f) {
		if (v2 & Matrices::R) {

			TimeEvent time(std::string("Update vector " + mNames(v1))); time.start();
			sumVectors(physics->instance()->f, physics->instance()->f, physics->instance()->R, alpha, beta);
			time.endWithBarrier(); _timeStatistics->addEvent(time);

			storeData(step, physics->instance()->f, "f", "right-hand side (f - R)");
			return;
		}
	}

	if (v1 & Matrices::B1c) {
		if (v2 & Matrices::primar) {

			TimeEvent time(std::string("Update vector " + mNames(v1))); time.start();
			#pragma omp parallel for
			for (size_t d = 0; d < physics->instance()->domains; d++) {
				for (size_t j = 0; j < physics->instance()->B1[d].J_col_indices.size(); j++) {
					if (physics->instance()->B1[d].I_row_indices[j] > (eslocal)physics->instance()->block[Instance::CONSTRAINT::DIRICHLET]) {
						break;
					}
					physics->instance()->B1c[d][j] = alpha * physics->instance()->B1c[d][j] + beta * physics->instance()->primalSolution[d][physics->instance()->B1[d].J_col_indices[j] - 1];
				}
			}
			time.endWithBarrier(); _timeStatistics->addEvent(time);

			storeData(step, physics->instance()->B1c, "B1c", "B1c (B1c - primar)");
			return;
		}
	}

	ESINFO(GLOBAL_ERROR) << "Implement updating vector " << mNames(v1);
}

void Solver::updateVector(const Step &step, Matrices v1, const std::vector<std::vector<double> > &v2, double alpha, double beta)
{
	ESINFO(PROGRESS2) << "Update vector: " << mNames(v1);
	if (v1 & Matrices::primar) {

		TimeEvent time(std::string("Update vector " + mNames(v1))); time.start();
		sumVectors(physics->instance()->primalSolution, physics->instance()->primalSolution, v2, alpha, beta);
		time.endWithBarrier(); _timeStatistics->addEvent(time);
	}
}

void Solver::regularizeMatrices(const Step &step, Matrices matrices)
{
	if (matrices & ~(Matrices::K)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid matrices passed to Solver::regularizeMatrices.";
	}

	ESINFO(PROGRESS2) << physics->name() << " regularizes matrices: " << mNames(matrices);

	TimeEvent time("Regularization of matrices " + mNames(matrices) + " by " + physics->name()); time.start();
	physics->makeStiffnessMatricesRegular(linearSolver->configuration.regularization);
	time.endWithBarrier(); _timeStatistics->addEvent(time);

	if (matrices & Matrices::K) {
		storeData(step, physics->instance()->N1, "N1", "N1");
		storeData(step, physics->instance()->N2, "N2", "N2");
		storeData(step, physics->instance()->RegMat, "RegMat", "RegMat");
	}
}

void Solver::composeGluing(const Step &step, Matrices matrices)
{
	if (matrices & ~(Matrices::B0 | Matrices::B1)) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid matrices passed to Solver::composeGluing.";
	}
	if (linearSolver->configuration.method == ESPRESO_METHOD::TOTAL_FETI) {
		matrices &= ~(Matrices::B0);
	}
	if (!matrices) {
		return;
	}

	ESINFO(PROGRESS2) << "Compose gluing matrices " << mNames(matrices) << "for " << physics->name();

	TimeEvent time("Composition of gluing matrices" + mNames(matrices) + "for " + physics->name()); time.start();
	if (matrices & Matrices::B0) {
		// TODO: create update method
		physics->instance()->B0.clear();
		physics->instance()->B0.resize(physics->instance()->domains);
		for (size_t d = 0; d < physics->instance()->domains; d++) {
			physics->instance()->B0[d].type = 'G';
			physics->instance()->B0subdomainsMap[d].clear();
		}
		switch (linearSolver->configuration.B0_type) {
		case B0_TYPE::CORNERS:
			physics->assembleB0FromCorners(step);
			break;
		case B0_TYPE::KERNELS:
			physics->assembleB0FromKernels(step);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown type of B0";
		}
	}
	if (matrices & Matrices::B1) {
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
	}
	time.endWithBarrier(); _timeStatistics->addEvent(time);

	if (matrices & Matrices::B0) {
		storeData(step, physics->instance()->B0, "B0", "B0");
	}

	if (matrices & Matrices::B1) {
		storeData(step, physics->instance()->B1, "B1", "B1");
		storeData(step, physics->instance()->B1c, "B1c", "B1c");
		storeData(step, physics->instance()->B1duplicity, "B1duplicity", "B1duplicity");
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

void Solver::initLinearSolver()
{
	ESINFO(PROGRESS2) << "Initialization of linear solver";

	TimeEvent timeSolver("Initialize linear solver"); timeSolver.startWithBarrier();
	linearSolver->init();
	timeSolver.end(); _timeStatistics->addEvent(timeSolver);
}

void Solver::updateLinearSolver(Matrices matrices)
{
	ESINFO(PROGRESS2) << "Updating of linear solver";

	TimeEvent timeSolver("Update linear solver"); timeSolver.startWithBarrier();
	linearSolver->update(matrices);
	timeSolver.end(); _timeStatistics->addEvent(timeSolver);
}

void Solver::runLinearSolver()
{
	ESINFO(PROGRESS2) << "Solve system";

	TimeEvent timeSolve("Linear Solver - runtime"); timeSolve.start();
	linearSolver->run();
	timeSolve.endWithBarrier(); _timeStatistics->addEvent(timeSolve);
}

void Solver::finalizeLinearSolver()
{
	ESINFO(PROGRESS2) << "Finalize " << _name << " solver";

	linearSolver->finilize();

	_store->finalize();

	_timeStatistics->totalTime.endWithBarrier();
	_timeStatistics->printStatsMPI();
}

void Solver::lineSearch(const std::vector<std::vector<double> > &U, std::vector<std::vector<double> > &deltaU, std::vector<std::vector<double> > &F_ext, Physics *physics, const Step &step)
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
		MPI_Allreduce(&cmul, &gmul, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return gmul;
	};

	double a = 0, b = 1, alpha = 1;
	double fa = 0, fb = 0, fx = 0, faStart = 0;

	std::vector<std::vector<double> > solution = deltaU;
	std::vector<std::vector<double> > F_ext_r = F_ext;

	for (size_t i = 0; i < 6; i++) {
		sumVectors(solution, U, deltaU, 1, alpha);

		solution.swap(physics->instance()->primalSolution);
		physics->assembleMatrix(step, Matrices::R);
		solution.swap(physics->instance()->primalSolution);

		if (i == 0) {
			faStart = multiply(deltaU, physics->instance()->f);
			sumVectors(F_ext_r, F_ext, physics->instance()->R, 1, -1);
			fb = multiply(deltaU, F_ext_r);
			if ((faStart < 0 && fb < 0) || (faStart >= 0 && fb >= 0)) {
				return;
			}
			fa = faStart;
		} else {
			sumVectors(F_ext_r, F_ext, physics->instance()->R, 1, -1);
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

	sumVectors(solution, U, deltaU, 0, alpha);
	solution.swap(deltaU);
}

void Solver::sumVectors(std::vector<std::vector<double> > &result, const std::vector<std::vector<double> > &a, const std::vector<std::vector<double> > &b, double alpha, double beta)
{
	#pragma omp parallel for
	for (size_t d = 0; d < a.size(); d++) {
		for (size_t i = 0; i < a[d].size(); i++) {
			result[d][i] = alpha * a[d][i] + beta * b[d][i];
		}
	}
}

