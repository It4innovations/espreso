
#include "apidataprovider.h"

#include "esinfo/eslog.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/meshinfo.h"

#include "config/reader/reader.h"

#include "mesh/mesh.h"
#include "output/resultstore.h"

#include "physics/system/builder/timebuilder.h"
#include "physics/kernels/heattransfer2d.kernel.h"
#include "physics/kernels/heattransfer3d.kernel.h"
#include "physics/kernels/structuralmechanics2d.kernel.h"
#include "physics/kernels/structuralmechanics3d.elasticity.kernel.h"
#include "physics/kernels/solverdataprovider/provider.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

APIDataProvider::APIDataProvider()
: kernel(NULL)
{
	info::system::setSignals();
	eslog::startln("ESPRESO: STARTED", "ESPRESO");

	ResultStore::createAsynchronizedStore();
}

APIDataProvider::~APIDataProvider()
{
	ResultStore::destroyAsynchronizedStore();
	eslog::endln("ESPRESO: FINISHED");
}

int APIDataProvider::nodesSize()
{
	return info::mesh->nodes->size;
}

int APIDataProvider::matrixType()
{
	return (int)kernel->solverDataProvider->feti->getMatrixType(0);
}

int APIDataProvider::DOFs()
{
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D: return 1; break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D: return 1; break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D: return 2; break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D: return 3; break;
	default: return 0;
	}
}

void APIDataProvider::prepare(int* argc, char ***argv)
{
	ECFReader::read(*info::ecf->ecfdescription, argc, argv, info::ecf->default_args, info::ecf->variables);
	info::ecf->input.decomposition.domains = 1;
	if (ResultStore::isComputeNode()) {
		Mesh::load();
		info::mesh->preprocess();
		info::mesh->printMeshStatistics();

		info::mesh->store->updateMesh();
		eslog::checkpointln("ESPRESO: MESH PREPARED");
	}

	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		kernel = new HeatTransfer2DKernel(NULL, info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d, info::ecf->heat_transfer_2d.load_steps_settings.at(step::loadstep + 1));
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		kernel = new HeatTransfer3DKernel(NULL, info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d, info::ecf->heat_transfer_3d.load_steps_settings.at(step::loadstep + 1));
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
		kernel = new StructuralMechanics2DKernel((StructuralMechanics2DKernel*)NULL, info::ecf->structural_mechanics_2d, info::ecf->structural_mechanics_2d, info::ecf->structural_mechanics_2d.load_steps_settings.at(step::loadstep + 1));
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		kernel = new StructuralMechanics3DKernel((StructuralMechanics3DKernel*)NULL, info::ecf->structural_mechanics_3d, info::ecf->structural_mechanics_3d, info::ecf->structural_mechanics_3d.load_steps_settings.at(step::loadstep + 1));
		break;
	default:
		eslog::globalerror("Physical solver: not implemented physical solver.\n");
	}
	info::mesh->store->updateMonitors();

	rhs.resize(DOFs() * nodesSize());
}

void APIDataProvider::fillMatrix(std::function<void(FETI4IInt, FETI4IInt, FETI4IInt*, FETI4IReal*)> add)
{
	FETI4IInt type, size, *nodes;
	TimeBuilder builder;
	Kernel::InstanceFiller filler(kernel->solutions.size());

	std::fill(rhs.begin(), rhs.end(), 0);
	builder.matrices = Builder::Request::K | Builder::Request::f;
	filler.begin = 0;
	filler.end = 1;
	filler.insert = [&] () {
		for (auto n = 0; n < size; ++n) {
			for (int d = 0; d < DOFs(); ++d) {
				rhs[DOFs() * nodes[n] + d] += filler.Fe[0][size * d + n];
			}
		}
		add(type, size, nodes, filler.Ke.vals);
	};
	auto enodes = info::mesh->elements->procNodes->begin();
	for (auto range = info::mesh->elements->eintervals.begin(); range != info::mesh->elements->eintervals.end(); ++range) {
		for (esint e = range->begin; e < range->end; ++e, filler.begin = filler.end++, ++enodes) {
			type = (int)Mesh::edata[range->code].type;
			size = Mesh::edata[range->code].nodes;
			nodes = enodes->data();
			kernel->processElements(builder, filler);
		}
	}
}

void APIDataProvider::fillDirichlet(std::vector<FETI4IInt> &indices, std::vector<FETI4IReal> &values)
{
	std::vector<std::pair<esint, esint> > dindices;
	kernel->solverDataProvider->general->dirichletIndices(dindices);
	std::vector<double > dvalues(dindices.size());
	kernel->solverDataProvider->general->dirichletValues(dvalues);

	std::vector<esint> permutation(dvalues.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		return dindices[i] < dindices[j];
	});

	indices.clear(); values.clear();
	for (size_t i = 0, j = 0; i < permutation.size(); i = j) {
		double sum = 0;
		while (j < permutation.size() && dindices[permutation[i]] == dindices[permutation[j]]) {
			sum += dvalues[permutation[j++]];
		}
		sum /= j - i;
		indices.push_back(dindices[permutation[i]].first * DOFs() + dindices[permutation[i]].second);
		values.push_back(sum);
	}
}

void APIDataProvider::fillL2G(std::vector<FETI4IInt> &l2g)
{
	l2g.assign(info::mesh->nodes->IDs->datatarray().begin(), info::mesh->nodes->IDs->datatarray().end());
}

void APIDataProvider::fillNeighbors(std::vector<FETI4IMPIInt> &neighbors)
{
	neighbors = info::mesh->neighbors;
}

void APIDataProvider::fillRHS(std::vector<FETI4IReal> &rhs)
{
	rhs = this->rhs;
}

void APIDataProvider::storeSolution(std::vector<FETI4IReal> &solution)
{
	memcpy(kernel->solutions.back().vals, solution.data(), sizeof(double) * solution.size());
	info::mesh->store->updateSolution();
}
