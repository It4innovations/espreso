
#include "mesio.h"

#include <vector>
#include <fstream>
#include <iomanip>

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

	MESIO mesio;
	MESIOInit(MPI_COMM_WORLD, 2);
	MESIOLoad(&mesio, MESIO_ENSIGHT, "benchmarks/input/cube/ensight/cube.case", MESIO_PARMETIS, 4);

	std::ofstream os("cube." + std::to_string(mpirank) + ".txt");
	{ // NODES
		MESIOInt nhalo, offset, size, totalSize;
		MESIOInt *ids, *position, *rankDist, *rankData;
		MESIOReal *coordinates;

		MESIONodes(mesio, &nhalo, &offset, &size, &totalSize, &ids, &position, &coordinates);
		MESIONodesRanks(mesio, &rankDist, &rankData);

		os << "NODES\n";
		os << "nhalo, offset, size / total size: " << nhalo << ", " << offset << ", " << size << " / " << totalSize << "\n";
		for (int n = 0; n < nhalo + size; ++n) {
			if (n == nhalo) {
				os << " -- my nodes -- \n";
			}
			os
			<< "N, p=" << std::setw(4) << position[n] << ", id=" << std::setw(4) << ids[n] << ", coo="
			<< std::setprecision(2) << std::scientific << coordinates[3 * n + 0] << ";"
			<< std::setprecision(2) << std::scientific << coordinates[3 * n + 1] << ";"
			<< std::setprecision(2) << std::scientific << coordinates[3 * n + 0] << ", ranks=";
			for (int r = 0, offset = rankDist[n]; r < mpisize; ++r) {
				if (offset < rankDist[n + 1] && rankData[offset] == r) {
					os << std::setw(2) << r;
					++offset;
				} else {
					os << "  ";
				}
			}
			os << "\n";
		}
	}

	{ // ELEMENTS
		MESIOInt offset, size, totalSize;
		MESIOInt *type, *enodesDist, *enodesData, *domains, *neighDist, *neighData;

		MESIOElements(mesio, &offset, &size, &totalSize, &type, &enodesDist, &enodesData);
		MESIOElementsDomains(mesio, &domains);
		MESIOElementsNeighbors(mesio, &neighDist, &neighData);

		os << "ELEMENTS\n";
		os << "offset, size / total size: " << offset << ", " << size << " / " << std::setw(4) << totalSize << "\n";
		for (int t = 0; t < MESIOElementType::SIZE; ++t) {
			MESIOElementsCounters(mesio, t, &offset, &totalSize);
			os << "type[" << std::setw(4) << t << "]: offset / total size: " << std::setw(4) << offset << " / " << std::setw(4) << totalSize << "\n";
		}
		for (int e = 0; e < size; ++e) {
			os << "E, offset=" << std::setw(4) << offset + e << ", domain=" << domains[e] << ", nodes[offset]=";
			for (int n = enodesDist[e]; n < enodesDist[e + 1]; ++n) {
				os << " " << std::setw(4) << enodesData[n];
			}
			os << "\n";
		}
		os << "DUAL GRAPH [-1 for external face]\n";
		for (int e = 0; e < size; ++e) {
			os << "DUAL element=" << std::setw(4) << offset + e << ", neighbors=";
			for (int n = neighDist[e]; n < neighDist[e + 1]; ++n) {
				os << " " << std::setw(4) << neighData[n];
			}
			os << "\n";
		}
		os << "\n";
	}

	{ // REGIONS OF ELEMENTS
		int regions = MESIOElementsRegions(mesio);
		for (int r = 0; r < regions; ++r) {
			{ // ELEMENTS
				const char* name;
				MESIOInt offset, size, totalSize, *elements;

				MESIOElementsRegion(mesio, r, &name, &size, &elements);
				os << "REGION OF ELEMENTS: " << name << "\n";
				for (int t = 0; t < MESIOElementType::SIZE; ++t) {
					MESIOElementsRegionCounters(mesio, r, t, &offset, &totalSize);
					os << "type[" << std::setw(4) << t << "]: offset / total size: " << std::setw(4) << offset << " / " << std::setw(4) << totalSize << "\n";
				}
				for (int n = 0; n < size; ++n) {
					os << "E, offset=" << std::setw(4) << elements[n] << "\n";
				}
			}
			{ // NODES
				MESIOInt nhalo, offset, size, totalSize, *nodes, *position;
				MESIOElementsRegionNodes(mesio, r, &nhalo, &offset, &size, &totalSize, &nodes, &position);

				os << "nhalo, offset, size / total size: " << nhalo << ", " << offset << ", " << size << " / " << totalSize << "\n";
				for (int n = 0; n < nhalo + size; ++n) {
					os << "N, p=" << std::setw(4) << position[n] << ", offset=" << std::setw(4) << nodes[n] << "\n";
				}
			}
		}
	}

	{ // BOUNDARY REGIONS
		int regions = MESIOBoundaryRegions(mesio);
		for (int r = 0; r < regions; ++r) {
			{ // ELEMENTS
				const char* name;
				MESIOInt dimension, offset, size, totalSize, *type, *parent, *enodesDist, *enodesData;

				MESIOBoundaryRegion(mesio, r, &name, &dimension, &size, &type, &parent, &enodesDist, &enodesData);
				os << "BOUNDARY REGION: " << name << "\n";
				if (dimension) {
					for (int t = 0; t < MESIOElementType::SIZE; ++t) {
						MESIOBoundaryRegionCounters(mesio, r, t, &offset, &totalSize);
						os << "type[" << std::setw(4) << t << "]: offset / total size: " << std::setw(4) << offset << " / " << std::setw(4) << totalSize << "\n";
					}
					for (int e = 0; e < size; ++e) {
						os << "E, offset=" << std::setw(4) << offset + e << ", parent=" << parent[e] << ", nodes[offset]=";
						for (int n = enodesDist[e]; n < enodesDist[e + 1]; ++n) {
							os << " " << std::setw(4) << enodesData[n];
						}
						os << "\n";
					}
				}
			}
			{ // NODES
				MESIOInt nhalo, offset, size, totalSize, *nodes, *position;
				MESIOBoundaryRegionNodes(mesio, r, &nhalo, &offset, &size, &totalSize, &nodes, &position);

				os << "nhalo, offset, size / total size: " << nhalo << ", " << offset << ", " << size << " / " << totalSize << "\n";
				for (int n = 0; n < nhalo + size; ++n) {
					os << "N, p=" << std::setw(4) << position[n] << ", offset=" << std::setw(4) << nodes[n] << "\n";
				}
			}
		}
	}

	MESIOFinalize();
	MPI_Finalize();
	return 0;
}




