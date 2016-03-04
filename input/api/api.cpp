
#include "api.h"

using namespace esinput;


void API::points(mesh::Coordinates &coordinates)
{
	mesh::Point p;

	eslocal max = 0;
	for (eslocal e = 0; e < eIndices.size(); e++) {
		max = std::max(max, *std::max_element(eIndices[e].begin(), eIndices[e].end()));
	}
	max /= DOFs;
	coordinates.reserve(max + 1);

	for (size_t i = 0; i <= max; i++) {
		coordinates.add(p, i, (ids[i * DOFs] - 1) / DOFs);
	}
}

void API::elements(std::vector<mesh::Element*> &elements)
{
	elements.reserve(eIndices.size());
	eslocal indices[20], params[6];

	for (eslocal e = 0; e < eIndices.size(); e++) {
		for (eslocal i = 0; i < eIndices[e].size(); i += DOFs) {
			indices[i / DOFs] = eIndices[e][i] / DOFs;
		}
		switch(eIndices[e].size() / DOFs) {
		case Tetrahedron4NodesCount:
			elements.push_back(new mesh::Tetrahedron4(indices, eIndices[e].size() / DOFs, params));
			break;
		case Tetrahedron10NodesCount:
			elements.push_back(new mesh::Tetrahedron10(indices, eIndices[e].size() / DOFs, params));
			break;
		case Pyramid5NodesCount:
			elements.push_back(new mesh::Pyramid5(indices, eIndices[e].size() / DOFs, params));
			break;
		case Pyramid13NodesCount:
			elements.push_back(new mesh::Pyramid13(indices, eIndices[e].size() / DOFs, params));
			break;
		case Prisma6NodesCount:
			elements.push_back(new mesh::Prisma6(indices, eIndices[e].size() / DOFs, params));
			break;
		case Prisma15NodesCount:
			elements.push_back(new mesh::Prisma15(indices, eIndices[e].size() / DOFs, params));
			break;
		case Hexahedron8NodesCount:
			elements.push_back(new mesh::Hexahedron8(indices, eIndices[e].size() / DOFs, params));
			break;
		case Hexahedron20NodesCount:
			elements.push_back(new mesh::Hexahedron20(indices, eIndices[e].size() / DOFs, params));
			break;
		default:
			ESLOG(eslog::ERROR) << "Unknown element with " << eIndices[e].size() / DOFs << " indices.";
		}
	}
}

void API::clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries)
{
	std::vector<esglobal> sBuffer;
	std::vector<std::vector<esglobal> > rBuffer(neighbours.size());
	std::vector<MPI_Request> req(2 * neighbours.size());
	std::vector<size_t> sizes(neighbours.size());

	size /= DOFs;

	for (size_t n = 0; n < neighbours.size(); n++) {
		MPI_Isend(&size,            sizeof(size_t), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n);
		MPI_Irecv(sizes.data() + n, sizeof(size_t), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + 2 * n + 1);
	}
	MPI_Waitall(2 * neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	sBuffer.reserve(size);
	for (size_t i = 0; i < size; i++) {
		sBuffer.push_back(ids[i * DOFs] / DOFs);
	}
	std::sort(sBuffer.begin(), sBuffer.end());

	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] != esconfig::MPIrank) {
			rBuffer[n].resize(sizes[n]);
		}
	}

	size_t rCounter = 0;
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] != esconfig::MPIrank) {
			MPI_Isend(sBuffer.data(),    sizes[n] * sizeof(esglobal), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + rCounter++);
			MPI_Irecv(rBuffer[n].data(), sizes[n] * sizeof(esglobal), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + rCounter++);
		}
	}
	MPI_Waitall(rCounter, req.data(), MPI_STATUSES_IGNORE);

	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, size);

	boundaries.resize(size);
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t n = 0; n < neighbours.size(); n++) {
				if (neighbours[n] == esconfig::MPIrank) {
					boundaries[i].push_back(esconfig::MPIrank);
				} else {
					auto it = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), ids[i * DOFs] / DOFs);
					if (it != rBuffer[n].end() && *it == ids[i * DOFs] / DOFs) {
						boundaries[i].push_back(neighbours[n]);
					}
				}
			}

		}
	}
}
