
#include "api.h"

using namespace espreso::input;


void API::points(Coordinates &coordinates)
{
	coordinates.reserveIndices(_size);
	for (size_t i = 0; i < _size; i++) {
		coordinates.add(i, _ids[i]);
	}
}

void API::elements(std::vector<Element*> &elements)
{
	elements.reserve(_eNodes.size());
	eslocal params[Element::PARAMS_SIZE];

	for (eslocal e = 0; e < _eNodes.size(); e++) {
		switch (_eType[e]) {
		case 0:
			ESTEST(MANDATORY) << "Point has to has only one index" << (_eNodes[e].size() != 1 ? TEST_FAILED : TEST_PASSED);
			elements.push_back(new UnknownPoint(_eNodes[e][0]));
			break;
		case 1:
			elements.push_back(new UnknownLine(_eNodes[e]));
			break;
		case 2:
			elements.push_back(new UnknownPlane(_eNodes[e]));
			break;
		case 3:
			elements.push_back(new UnknownVolume(_eNodes[e]));
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown element type " << _eType[e];
		}
	}
}

void API::clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours)
{
	_mesh.DOFs().reserve(_size);
	for (size_t i = 0; i < _size; i++) {
		_mesh.DOFs().push_back(new DOF(i));
	}
	auto it = std::find(_neighbours.begin(), _neighbours.end(), config::env::MPIrank);
	if (it != _neighbours.end() && *it == config::env::MPIrank) {
		_neighbours.erase(it);
	}
	std::vector<eslocal> sBuffer;
	std::vector<std::vector<eslocal> > rBuffer(_neighbours.size());
	std::vector<MPI_Request> req(2 * _neighbours.size());
	std::vector<size_t> sizes(_neighbours.size());

	for (size_t n = 0; n < _neighbours.size(); n++) {
		MPI_Isend(&_size,           sizeof(size_t), MPI_BYTE, _neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n);
		MPI_Irecv(sizes.data() + n, sizeof(size_t), MPI_BYTE, _neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n + 1);
	}
	MPI_Waitall(2 * _neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	sBuffer.reserve(_size);
	for (size_t i = 0; i < _size; i++) {
		sBuffer.push_back(_ids[i]);
	}
	std::sort(sBuffer.begin(), sBuffer.end());

	for (size_t n = 0; n < _neighbours.size(); n++) {
		rBuffer[n].resize(sizes[n]);
	}

	for (size_t n = 0; n < _neighbours.size(); n++) {
		MPI_Isend(sBuffer.data(),       _size * sizeof(eslocal), MPI_BYTE, _neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n);
		MPI_Irecv(rBuffer[n].data(), sizes[n] * sizeof(eslocal), MPI_BYTE, _neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n + 1);
	}
	MPI_Waitall(2 * _neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _size);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t n = 0; n < _neighbours.size(); n++) {
				if (n == pushMyRank) {
					nodes[i]->clusters().push_back(config::env::MPIrank);
				}
//				auto it = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), _ids[i * _DOFs] / _DOFs);
//				if (it != rBuffer[n].end() && *it == _ids[i * _DOFs] / _DOFs) {
//					nodes[i]->clusters().push_back(_neighbours[n]);
//					realNeighbour[t].insert(_neighbours[n]);
//				}
			}
			if (_neighbours.size() == pushMyRank) {
				nodes[i]->clusters().push_back(config::env::MPIrank);
			}

		}
	}

	for (size_t t = 1; t < threads; t++) {
		realNeighbour[0].insert(realNeighbour[t].begin(), realNeighbour[t].end());
	}

	neighbours = std::vector<int>(realNeighbour[0].begin(), realNeighbour[0].end());
}
