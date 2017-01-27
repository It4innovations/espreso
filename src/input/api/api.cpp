
#include "api.h"
#include "../../config/environment.h"
#include "../../mesh/settings/evaluator.h"

#include "../../mesh/elements/volume/unknownvolume.h"
#include "../../mesh/elements/plane/unknownplane.h"
#include "../../mesh/elements/line/unknownline.h"
#include "../../mesh/elements/point/unknownpoint.h"
#include "../../mesh/elements/point/dof.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"

#include "../../basis/utilities/utils.h"

using namespace espreso::input;


void API::points(const std::vector<std::vector<eslocal> > &eNodes, size_t DOFsSize)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, eNodes.size());

	std::vector<eslocal> tMax(threads);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		eslocal max = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			max = std::max(max, *Esutils::max_element(eNodes[e]));
		}
		tMax[t] = max;
	}

	_mesh._coordinates->resize(*std::max_element(tMax.begin(), tMax.end()) + 1);

	_mesh._DOFs.reserve(DOFsSize);
	for (size_t d = 0; d < DOFsSize; d++) {
		_mesh._DOFs.push_back(new DOF(d));
	}

	_mesh.fillNodesFromCoordinates();
}

void API::elements(const std::vector<eslocal> &eType, std::vector<std::vector<eslocal> > &eNodes, std::vector<std::vector<eslocal> > &eDOFs, std::vector<std::vector<double> > &eMatrices)
{
	_mesh._elements.reserve(eNodes.size());

	for (size_t e = 0; e < eNodes.size(); e++) {
		switch (eType[e]) {
		case 0:
			ESTEST(MANDATORY) << "Point has to has only one index" << (eNodes[e].size() != 1 ? TEST_FAILED : TEST_PASSED);
			_mesh._elements.push_back(new UnknownPoint(eNodes[e][0]));
			break;
		case 1:
			_mesh._elements.push_back(new UnknownLine(_mesh.nodes(), eNodes[e], eDOFs[e], eMatrices[e]));
			break;
		case 2:
			_mesh._elements.push_back(new UnknownPlane(_mesh.nodes(), eNodes[e], eDOFs[e], eMatrices[e]));
			break;
		case 3:
			_mesh._elements.push_back(new UnknownVolume(_mesh.nodes(), eNodes[e], eDOFs[e], eMatrices[e]));
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown element type " << eType[e];
		}
	}

	_mesh.fillParentElementsToNodes();
	_mesh.fillParentElementsToDOFs(eDOFs);

	_mesh.partitiate(_configuration.domains);
}

void API::dirichlet(size_t dirichletSize, eslocal *dirichletIndices, double *dirichletValues)
{
	_mesh._evaluators.push_back(new ArrayEvaluator("dirichletAPI", dirichletSize, dirichletIndices, dirichletValues, _offset, Property::UNKNOWN));
	_mesh._regions.push_back(new Region());

	_mesh._regions[0]->settings.resize(1);
	_mesh._regions[0]->settings[0][Property::UNKNOWN].push_back(_mesh._evaluators.back());
	_mesh._regions[0]->elements().resize(dirichletSize);

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, dirichletSize);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			_mesh._DOFs[dirichletIndices[i] - _offset]->regions().push_back(_mesh._regions[0]);
			_mesh._regions[0]->elements()[i] = _mesh._DOFs[dirichletIndices[i] - _offset];
		}
	}
}

void API::clusterBoundaries(std::vector<int> &neighbours, size_t size, const eslocal *l2g)
{
	auto it = std::find(neighbours.begin(), neighbours.end(), environment->MPIrank);
	if (it != neighbours.end() && *it == environment->MPIrank) {
		neighbours.erase(it);
	}

	std::vector<std::vector<eslocal> > rBuffer(neighbours.size());
	std::vector<MPI_Request> req(2 * neighbours.size());
	std::vector<size_t> sizes(neighbours.size());

	for (size_t n = 0; n < neighbours.size(); n++) {
		MPI_Isend(&size           , sizeof(size_t), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n);
		MPI_Irecv(sizes.data() + n, sizeof(size_t), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n + 1);
	}
	MPI_Waitall(2 * neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	std::vector<eslocal> sBuffer(l2g, l2g + size);
	std::sort(sBuffer.begin(), sBuffer.end());

	for (size_t n = 0; n < neighbours.size(); n++) {
		rBuffer[n].resize(sizes[n]);
	}

	for (size_t n = 0; n < neighbours.size(); n++) {
		MPI_Isend(sBuffer.data(),        size * sizeof(eslocal), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n);
		MPI_Irecv(rBuffer[n].data(), sizes[n] * sizeof(eslocal), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n + 1);
	}
	MPI_Waitall(2 * neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, size);

	size_t pushMyRank = std::lower_bound(neighbours.begin(), neighbours.end(), environment->MPIrank) - neighbours.begin();
	std::vector<std::vector<G2L> > g2l(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t n = 0; n < neighbours.size(); n++) {
				if (n == pushMyRank) {
					_mesh._DOFs[i]->clusters().push_back(environment->MPIrank);
				}
				auto it = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), l2g[i]);
				if (it != rBuffer[n].end() && *it == l2g[i]) {
					_mesh._DOFs[i]->clusters().push_back(neighbours[n]);
				}
			}
			if (neighbours.size() == pushMyRank) {
				_mesh._DOFs[i]->clusters().push_back(environment->MPIrank);
			}
			if (_mesh._DOFs[i]->clusters().size() > 1) {
				g2l[t].push_back(G2L(l2g[i], i));
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		_mesh._g2l->insert(_mesh._g2l->end(), g2l[t].begin(), g2l[t].end());
	}
	std::sort(_mesh._g2l->begin(), _mesh._g2l->end());

	_mesh._neighbours = neighbours;
}
