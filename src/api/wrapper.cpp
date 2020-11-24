
#include "feti4i.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/timelogger.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/sysutils.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "config/configuration.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/preprocessing/meshpreprocessing.h"

#include "physics/system/fetisystem.h"
#include "physics/system/builder/timebuilder.h"
#include "physics/composer/feti/nodes.uniform.api.composer.h"
#include "feti/generic/FETISystemSolver.h"

#include "omp.h"

#include <numeric>
#include <vector>
#include <unordered_map>

struct FETI4IStructMatrix {
	FETI4IStructMatrix(esint size, esint *l2g, esint indexing, esint type, esint dofs, esint reorder)
	: size(size), indexing(indexing), type(type), dofs(dofs), reorder(reorder),
	  l2g(l2g, l2g + size)
	{
		int threads = espreso::info::env::OMP_NUM_THREADS;
		etype.resize(threads);
		nodes.resize(threads);
		stiffness.resize(threads);

		ndist.resize(threads);
		sdist.resize(threads);

		#pragma omp parallel for
		for (int t = 0; t < threads; ++t) {
			etype[t] = new std::vector<espreso::Element*>();
			nodes[t] = new std::vector<esint>();
			stiffness[t] = new std::vector<double>();

			ndist[t] = new std::vector<esint>();
			sdist[t] = new std::vector<esint>();
		}
		ndist[0]->push_back(0);
		sdist[0]->push_back(0);
	};

	espreso::Mesh mesh;

	esint size;
	esint indexing;
	esint type;
	esint dofs;
	esint reorder;

	std::vector<esint> l2g;

	std::vector<std::vector<espreso::Element*>* > etype;
	std::vector<std::vector<esint>* > nodes, ndist, sdist;
	std::vector<std::vector<double>* > stiffness;
};

struct FETI4IStructInstance {
	FETI4IStructMatrix *matrix;
	espreso::FETISolverData data;
	espreso::FETISystemSolver solver;

	FETI4IStructInstance(FETI4IStructMatrix *matrix)
	: matrix(matrix),
	  data(espreso::info::ecf->feti4ilibrary.solver),
	  solver(espreso::info::ecf->feti4ilibrary.solver, data)
	{

	}
};

struct FETI4IDataHolder {
	std::unordered_map<void*, FETI4IStructMatrix*> matrices;
	std::unordered_map<void*, FETI4IStructInstance*> instances;
};

FETI4IDataHolder FETI4IData;

using namespace espreso;

void error(const std::string &error)
{
	printf("error: %s\n", error.c_str());
}

void FETI4IInit(
		MPI_Comm 		comm,
		FETI4IInt		verbosity)
{
	info::env::set();
	info::mpi::init(comm);
	MPITools::init();
	eslog::init(new Logger<ProgressTerminalLogger>);
	if (utils::exists("espreso.ecf")) {
		ECF::init("espreso.ecf");
	} else {
		ECF::init();
	}
	for (int i = 0; i < eslog::logger->size; ++i) {
		eslog::logger->args[i]->verbosity = verbosity;
	}
	Mesh::init();

	eslog::startln("FETI4I: INITIALIZED", "FETI4I");
}

void FETI4IFinalize()
{
	Mesh::finish();
	ECF::finish();
	MPITools::finish();
	eslog::endln("FETI4I: FINISHED");
}

void FETI4ISetDefaultIntegerOptions(
		FETI4IInt*		options)
{
	const auto &settings = info::ecf->feti4ilibrary;
	options[FETI4I_SUBDOMAINS] = settings.domains;

	options[FETI4I_MAX_ITERATIONS] = settings.solver.max_iterations;
	options[FETI4I_FETI_METHOD] = static_cast<int>(settings.solver.method);
	options[FETI4I_PRECONDITIONER] = static_cast<int>(settings.solver.preconditioner);
	options[FETI4I_CGSOLVER] = static_cast<int>(settings.solver.iterative_solver);
}

void FETI4ISetDefaultRealOptions(
		FETI4IReal*		options)
{
	const auto &settings = info::ecf->feti4ilibrary;
	options[FETI4I_PRECISION] = settings.solver.precision;
}

void FETI4ICreateStiffnessMatrix(
		FETI4IMatrix	*matrix,
		FETI4IInt		nodesSize,
		FETI4IInt*		l2g,
		FETI4IInt		indexBase,
		FETI4IInt		type,
		FETI4IInt		dofsPerNode,
		FETI4IInt		dofsOrdering)
{
	*matrix = new FETI4IStructMatrix(nodesSize, l2g, indexBase, type, dofsPerNode, dofsOrdering);
	FETI4IData.matrices[*matrix] = *matrix;
}

void FETI4IAddElement(
		FETI4IMatrix	matrix,
		FETI4IInt		type,
		FETI4IInt		size,
		FETI4IInt*		nodes,
		FETI4IReal*		stiffness)
{
	int thread = omp_get_thread_num();
	switch (type) {
	case FETI4I_ETYPE_POINT:
		switch (size) {
		case 1: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::POINT1]); break;
		default: error("invalid element type / size combination");
		}
		break;
	case FETI4I_ETYPE_LINE:
		switch (size) {
		case 2: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::LINE2]); break;
		case 3: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::LINE3]); break;
		default: error("invalid element type / size combination");
		}
		break;
	case FETI4I_ETYPE_PLANE:
		switch (size) {
		case 3: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::TRIANGLE3]); break;
		case 4: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::SQUARE4]); break;
		case 6: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::TRIANGLE6]); break;
		case 8: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::SQUARE8]); break;
		default: error("invalid element type / size combination");
		}
		break;
	case FETI4I_ETYPE_VOLUME:
		switch (size) {
		case  4: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::TETRA4]); break;
		case  5: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::PYRAMID5]); break;
		case  6: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::PRISMA6]); break;
		case  8: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::HEXA8]); break;
		case 10: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::TETRA10]); break;
		case 13: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::PYRAMID13]); break;
		case 15: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::PRISMA15]); break;
		case 20: matrix->etype[thread]->push_back(&Mesh::edata[(int)Element::CODE::HEXA20]); break;
		default: error("invalid element type / size combination");
		}
		break;
	}

	if (matrix->reorder) {
		for (FETI4IInt row = 0; row < size; ++row) {
			matrix->nodes[thread]->push_back(nodes[row] - matrix->indexing);
			for (FETI4IInt rdof = 0; rdof < matrix->dofs; ++rdof) {
				for (FETI4IInt col = 0; col < size; ++col) {
					for (FETI4IInt cdof = 0; cdof < matrix->dofs; ++cdof) {
						FETI4IInt r = row + size * rdof;
						FETI4IInt c = col + size * cdof;
						matrix->stiffness[thread]->push_back(stiffness[r * (size * matrix->dofs) + c]);
					}
				}
			}
		}
	} else {
		for (FETI4IInt row = 0; row < size; ++row) {
			matrix->nodes[thread]->push_back(nodes[row] - matrix->indexing);
			for (FETI4IInt rdof = 0; rdof < matrix->dofs; ++rdof) {
				for (FETI4IInt col = 0; col < size; ++col) {
					for (FETI4IInt cdof = 0; cdof < matrix->dofs; ++cdof) {
						FETI4IInt r = ((row * matrix->dofs) + rdof);
						FETI4IInt c = ((col * matrix->dofs) + cdof);
						matrix->stiffness[thread]->push_back(stiffness[r * (size * matrix->dofs) + c]);
					}
				}
			}
		}
	}

	matrix->ndist[thread]->push_back(matrix->nodes[thread]->size());
	matrix->sdist[thread]->push_back(matrix->stiffness[thread]->size());
}

void FETI4ICreateInstance(
		FETI4IInstance	*instance,
		FETI4IMatrix	matrix,
		FETI4IMPIInt	neighbours_size,
		FETI4IMPIInt*	neighbours,
		FETI4IInt		dirichlet_size,
		FETI4IInt*		dirichlet_indices,
		FETI4IReal*		dirichlet_values,
		FETI4IInt*		iopts,
		FETI4IReal*		ropts)
{
	int threads = espreso::info::env::OMP_NUM_THREADS;

	auto &settings = info::ecf->feti4ilibrary;

	settings.domains = iopts[FETI4I_SUBDOMAINS];
	settings.solver.max_iterations = iopts[FETI4I_MAX_ITERATIONS];
	settings.solver.ecfdescription->getParameter(&settings.solver.method)->setValue(std::to_string(iopts[FETI4I_FETI_METHOD]));
	settings.solver.ecfdescription->getParameter(&settings.solver.preconditioner)->setValue(std::to_string(iopts[FETI4I_PRECONDITIONER]));
	settings.solver.ecfdescription->getParameter(&settings.solver.iterative_solver)->setValue(std::to_string(iopts[FETI4I_CGSOLVER]));
	settings.solver.precision = ropts[FETI4I_PRECISION];

	FETI4IInstance system = *instance = new FETI4IStructInstance(matrix);
	FETI4IData.instances[system] = system;

	std::vector<std::vector<Element*> > etype(threads);
	std::vector<std::vector<esint> > nodes(threads), ndist(threads), sdist(threads);
	std::vector<std::vector<double> > stiffness(threads);
	for (int t = 0; t < threads; ++t) {
		matrix->etype[t]->swap(etype[t]);
		matrix->nodes[t]->swap(nodes[t]);
		matrix->ndist[t]->swap(ndist[t]);
		matrix->sdist[t]->swap(sdist[t]);
		matrix->stiffness[t]->swap(stiffness[t]);
	}

	utils::threadDistributionToFullDistribution(ndist);
	utils::threadDistributionToFullDistribution(sdist);

	if (!std::is_sorted(matrix->l2g.begin(), matrix->l2g.end())) {
		std::vector<esint> permutation(matrix->l2g.size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
			return matrix->l2g[i] < matrix->l2g[j];
		});

		std::sort(matrix->l2g.begin(), matrix->l2g.end());

		std::vector<esint> backpermutation(permutation.size());
		std::iota(backpermutation.begin(), backpermutation.end(), 0);
		std::sort(backpermutation.begin(), backpermutation.end(), [&] (esint i, esint j) { return permutation[i] < permutation[j]; });

		for (size_t t = 0; t < nodes.size(); ++t) {
			for (size_t n = 0; n < nodes[t].size(); ++n) {
				nodes[t][n] = backpermutation[nodes[t][n]];
			}
		}
	}

	matrix->mesh.neighbors = std::vector<int>(neighbours, neighbours + neighbours_size);
	matrix->mesh.neighborsWithMe = matrix->mesh.neighbors;
	matrix->mesh.neighborsWithMe.push_back(info::mpi::rank);
	std::sort(matrix->mesh.neighborsWithMe.begin(), matrix->mesh.neighborsWithMe.end());
	matrix->mesh.nodes->size = matrix->l2g.size();
	matrix->mesh.nodes->IDs = new serializededata<esint, esint>(1, tarray<esint>(threads, matrix->l2g));
	matrix->mesh.nodes->distribution = matrix->mesh.nodes->IDs->datatarray().distribution();
	matrix->mesh.elements->epointers = new serializededata<esint, Element*>(1, etype);
	matrix->mesh.elements->procNodes = new serializededata<esint, esint>(ndist, nodes);
	matrix->mesh.elements->stiffness = new serializededata<esint, double>(sdist, stiffness);

	matrix->mesh.elements->offset = matrix->mesh.elements->epointers->datatarray().size();;
	matrix->mesh.elements->size = matrix->mesh.elements->epointers->datatarray().size();
	matrix->mesh.elements->totalSize = Communication::exscan(matrix->mesh.elements->offset);
	matrix->mesh.elements->distribution = matrix->mesh.elements->epointers->datatarray().distribution();
	matrix->mesh.elements->IDs = new serializededata<esint, esint>(1, tarray<esint>(matrix->mesh.elements->distribution, 1));
	std::iota(matrix->mesh.elements->IDs->datatarray().begin(), matrix->mesh.elements->IDs->datatarray().end(), matrix->mesh.elements->offset);

	Mesh *mesh = &matrix->mesh;
	std::swap(mesh, info::mesh);
	mesh::computeNodesDuplication();
	mesh::partitiate(settings.domains, false);
	if (settings.solver.method == FETIConfiguration::METHOD::HYBRID_FETI) {
		mesh::computeDomainDual();
	}

	NodesUniformAPIComposer composer(settings.solver, matrix->dofs);
	composer.fill(system->data);

	struct FETI4IIJV { FETI4IInt i, j; FETI4IReal v; };

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<esint> dofs;
		for (esint d = info::mesh->elements->domainDistribution[t]; d != info::mesh->elements->domainDistribution[t + 1]; ++d) {
			system->data.K[d].type = static_cast<MatrixType>(matrix->type);

			std::vector<FETI4IIJV> m;
			auto enodes = info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d];
			auto stiffness = info::mesh->elements->stiffness->begin() + info::mesh->elements->elementsDistribution[d];
			for (esint e = info::mesh->elements->elementsDistribution[d]; e < info::mesh->elements->elementsDistribution[d + 1]; ++e, ++enodes, ++stiffness) {
				dofs.clear();
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					auto dmap = composer.DOFMap()->begin() + (*n * matrix->dofs);
					for (int dof = 0; dof < matrix->dofs; ++dof, ++dmap) {
						for (auto di = dmap->begin(); di != dmap->end(); ++di) {
							if (di->domain == d + info::mesh->elements->firstDomain) {
								dofs.push_back(di->index);
							}
						}
					}
				}
				for (size_t i = 0; i < dofs.size(); ++i) {
					for (size_t j = 0; j < dofs.size(); ++j) {
						if (dofs[i] <= dofs[j] || matrix->type == 2) {
							m.push_back(FETI4IIJV{ dofs[i] + 1, dofs[j] + 1, stiffness->data()[i * dofs.size() + j]});
						}
					}
				}
			}

			std::sort(m.begin(), m.end(), [] (FETI4IIJV &v1, FETI4IIJV &v2) {
				if (v1.i == v2.i) {
					return v1.j < v2.j;
				}
				return v1.i < v2.i;
			});

			std::vector<esint> rows({ 1 }), cols = { m.front().j };
			std::vector<double> vals = { m.front().v };
			for (size_t i = 1, nonzeros = 0; i < m.size(); i++) {
				if (m[i].i != m[i - 1].i || m[i].j != m[i - 1].j) {
					++nonzeros;
					cols.push_back(m[i].j);
					vals.push_back(m[i].v);
					if (m[i - 1].i != m[i].i) {
						rows.push_back(nonzeros + 1);
					}
				} else {
					vals.back() += m[i].v;
				}
			}
			rows.push_back(cols.size() + 1);

			system->data.K[d].resize(rows.size() - 1, rows.size() - 1, cols.size());
			system->data.K[d].fillPattern(rows.size() - 1, rows.data(), cols.data());
			system->data.K[d].structureUpdated();
			system->data.K[d].fillValues(vals.size(), vals.data());
			system->data.f.resizeDomain(d, rows.size() - 1);
		}
	}


	system->data.x.shallowCopyStructure(&system->data.f);
	system->data.x.setDuplications(DataDecomposition::DUPLICATION::DUPLICATE);
	system->data.y.shallowCopyStructure(&system->data.f);
	system->data.y.setDuplications(DataDecomposition::DUPLICATION::DUPLICATE);

	system->data.BC.initVectors(1);
	system->data.BC.resize(info::mesh->nodes->size * matrix->dofs, dirichlet_size);
	system->data.BC.holder()->fillPattern(dirichlet_indices);
	system->data.BC.holder()->fillValues(dirichlet_values);
	system->data.Kdiag.shallowCopyStructure(system->data.f.at(0));

	TimeBuilder builder;
	builder.matrices = Builder::Request::KCM | Builder::Request::RBCf;
	system->data.buildB1();
	system->data.setDirichlet(&builder);

	settings.solver.regularization = FETIConfiguration::REGULARIZATION::ALGEBRAIC;
	settings.solver.B0_type = FETIConfiguration::B0_TYPE::KERNELS;
	system->solver.init();
	system->solver.update();

	std::swap(info::mesh, mesh);
}

void FETI4ISolve(
		FETI4IInstance	instance,
		FETI4IReal*		rhs,
		FETI4IReal*		solution)
{
	VectorDense f(instance->matrix->size * instance->matrix->dofs, rhs), result(instance->matrix->size * instance->matrix->dofs, solution);
	instance->data.f[0].fillData(&f);
	instance->solver.insertRHS(instance->data.f);

	TimeBuilder builder;
	builder.matrices = Builder::Request::KCM | Builder::Request::RBCf;
	instance->data.printData(&builder, ("api/" + std::to_string(info::mpi::rank)).c_str());

	instance->solver.solve();
	result.fillData(&instance->data.x[0]);
}

void FETI4IUpdateStiffnessMatrix(
		FETI4IInstance	instance,
		FETI4IMatrix	stiffnessMatrix)
{
	error("not imlemented method\n");
}

void FETI4IUpdateRhs(
		FETI4IInstance	instance,
		FETI4IInt		size,
		FETI4IReal*		values)
{
	error("not imlemented method\n");
}

void FETI4IUpdateDirichlet(
		FETI4IInstance	instance,
		FETI4IInt		size,
		FETI4IInt*		indices,
		FETI4IReal*		values)
{
	error("not imlemented method\n");
}

void FETI4IDestroy(
		void*			ptr)
{
	auto mit = FETI4IData.matrices.find(ptr);
	if (mit != FETI4IData.matrices.end()) {
		delete mit->second;
	}
	auto iit = FETI4IData.instances.find(ptr);
	if (iit != FETI4IData.instances.end()) {
		delete iit->second;
	}
}
