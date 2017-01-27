
#include "esdata.h"

#include "../../config/environment.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/elements/element.h"

using namespace espreso::store;

void Esdata::mesh(const Mesh &mesh, const std::string &path)
{
	Esdata(mesh, path);
}

Esdata::Esdata(const Mesh &mesh, const std::string &path)
: _mesh(mesh), _path(path)
{
	std::stringstream ss;
	ss << "mkdir -p " << _path;
	int out = system(ss.str().c_str());
	for (size_t p = 0; p < _mesh.parts(); p++) {
		std::stringstream ssDir;
		ssDir << ss.str() << "/" << p + _mesh.parts() * environment->MPIrank;
		out = system(ssDir.str().c_str());
	}

	if (out) {
		ESINFO(ERROR) << "Cannot create output directory";
	}

	coordinates(_mesh.coordinates());
	elements(_mesh);
	regions(_mesh);
	boundaries(_mesh);
}

void Esdata::coordinates(const Coordinates &coordinates)
{
	#pragma omp parallel for
	for  (size_t p = 0; p < coordinates.parts(); p++) {
		std::ofstream os;
		eslocal size;
		esglobal index;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/coordinates.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		size = coordinates.localSize(p);
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (size_t i = 0; i < coordinates.localSize(p); i++) {
			index = coordinates.globalIndex(i, p);
			os.write(reinterpret_cast<const char*>(&index), sizeof(esglobal));
			const Point &point = coordinates.get(i, p);
			os.write(reinterpret_cast<const char*>(&point), Point::size() * sizeof(double));
		}
		os.close();
	}
}

void Esdata::elements(const Mesh &mesh)
{
	#pragma omp parallel for
	for  (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		const std::vector<eslocal> &parts = mesh.getPartition();
		const std::vector<Element*> &elements = mesh.elements();
		eslocal size;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/elements.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		// elements
		size = parts[p + 1] - parts[p];
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			elements[e]->store(os, mesh.coordinates(), p);
		}

		os.close();
	}
}

//void Esdata::materials(const Mesh &mesh, const std::vector<Material*> &materials)
//{
//	#pragma omp parallel for
//	for (size_t p = 0; p < mesh.parts(); p++) {
//		std::ofstream os;
//		std::stringstream ss;
//		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/materials.dat";
//		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
//
//		eslocal size = materials.size();
//		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
//		for (size_t i = 0; i < materials.size(); i++) {
//			os << materials[i];
//		}
//		os.close();
//	}
//}

void Esdata::regions(const Mesh &mesh)
{
	#pragma omp parallel for
	for  (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/settings.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		eslocal size;
		// regions
		size = _mesh.regions().size();
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (size_t i = 0; i < _mesh.regions().size(); i++) {
			eslocal length = _mesh.regions()[i]->name.size();
			os.write(reinterpret_cast<const char *>(&length), sizeof(eslocal));
			os.write(_mesh.regions()[i]->name.c_str(), _mesh.regions()[i]->name.size());
		}
		os.close();
	}
	ESINFO(GLOBAL_ERROR) << "Implement store regions";
}

void Esdata::boundaries(const Mesh &mesh)
{
//	Boundaries boundaries = mesh.subdomainBoundaries();
//	const Boundaries &cBoundaries = mesh.clusterBoundaries();
//
//	size_t threads = environment->OMP_NUM_THREADS;
//	std::vector<size_t> distribution = Esutils::getDistribution(threads, boundaries.size());
//
//	std::vector<std::vector<std::vector<esglobal> > > sBuffer(threads, std::vector<std::vector<esglobal> > (mesh.neighbours().size()));
//	std::vector<std::vector<esglobal> > rBuffer(mesh.neighbours().size());
//
//	cilk_for (size_t t = 0; t < threads; t++) {
//		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
//			for (size_t j = 0; j < boundaries[i].size(); j++) {
//				boundaries[i][j] += environment->MPIrank * mesh.parts();
//			}
//			for (size_t n = 0; n < mesh.neighbours().size(); n++) {
//				if (std::binary_search(cBoundaries[i].begin(), cBoundaries[i].end(), mesh.neighbours()[n])) {
//					sBuffer[t][n].push_back(mesh.coordinates().globalIndex(i));
//					sBuffer[t][n].push_back(boundaries[i].size());
//					sBuffer[t][n].insert(sBuffer[t][n].end(), boundaries[i].begin(), boundaries[i].end());
//				}
//			}
//		}
//	}
//
//	cilk_for (size_t n = 0; n < mesh.neighbours().size(); n++) {
//		for (size_t t = 1; t < threads; t++) {
//			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
//		}
//	}
//
//	std::vector<MPI_Request> req(mesh.neighbours().size());
//
//	for (size_t n = 0; n < mesh.neighbours().size(); n++) {
//		MPI_Isend(sBuffer[0][n].data(), sizeof(esglobal) * sBuffer[0][n].size(), MPI_BYTE, mesh.neighbours()[n], 0, MPI_COMM_WORLD, req.data() + n);
//	}
//
//	auto n2i = [ & ] (size_t neighbour) {
//		return std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), neighbour) - mesh.neighbours().begin();
//	};
//
//	int flag, counter = 0;
//	MPI_Status status;
//	while (counter < mesh.neighbours().size()) {
//		MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
//		if (flag) {
//			int count;
//			MPI_Get_count(&status, MPI_BYTE, &count);
//			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(esglobal));
//			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//			counter++;
//		}
//	}
//
//	MPI_Waitall(mesh.neighbours().size(), req.data(), MPI_STATUSES_IGNORE);
//
//	for (size_t n = 0; n < mesh.neighbours().size(); n++) {
//		for (size_t i = 0; i < rBuffer[n].size(); i++) {
//			esglobal index = mesh.coordinates().clusterIndex(rBuffer[n][i]);
//			boundaries[index].insert(boundaries[index].end(), rBuffer[n].begin() + i + 2, rBuffer[n].begin() + i + 2 + rBuffer[n][i + 1]);
//			i += rBuffer[n][i + 1] + 1;
//		}
//	}
//
//	cilk_for (size_t t = 0; t < threads; t++) {
//		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
//			std::sort(boundaries[i].begin(), boundaries[i].end());
//		}
//	}
//
//	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
//		std::ofstream os;
//		std::stringstream ss;
//		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/boundaries.dat";
//		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
//
//		eslocal value, size;
//		esglobal index;
//		size = 0;
//		for (size_t i = 0; i < boundaries.size(); i++) {
//			if (std::binary_search(boundaries[i].begin(), boundaries[i].end(), p + _mesh.parts() * environment->MPIrank)) {
//				size++;
//			}
//		}
//		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
//
//		for (size_t i = 0; i < boundaries.size(); i++) {
//			if (std::binary_search(boundaries[i].begin(), boundaries[i].end(), p + _mesh.parts() * environment->MPIrank)) {
//				size = boundaries[i].size();
//				os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
//				for (auto it = boundaries[i].begin(); it != boundaries[i].end(); ++it) {
//					value = *it;
//					os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
//				}
//			}
//		}
//		os.close();
//	}

	#pragma omp parallel for
	for  (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/boundaries.dat";
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		eslocal domain, size;

		size = mesh.coordinates().localSize(p);
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));

		for (size_t i = 0; i < mesh.coordinates().localSize(p); i++) {
			const std::vector<eslocal> &domains = mesh.nodes()[mesh.coordinates().clusterIndex(i, p)]->domains();

			eslocal dSize = domains.size();
			os.write(reinterpret_cast<const char*>(&dSize), sizeof(eslocal));
			for (size_t d = 0; d < domains.size(); d++) {
				domain = domains[d];
				os.write(reinterpret_cast<const char*>(&domain), sizeof(eslocal));
			}
		}

		os.close();
	}

}

