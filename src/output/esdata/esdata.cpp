
#include "esdata.h"

using namespace espreso::output;

void Esdata::store(double shrinkSubdomain, double shrinkCluster)
{
	// shrink is ignored!!!

	std::stringstream ss;
	ss << "mkdir -p " << _path;
	system(ss.str().c_str());
	for (size_t p = 0; p < _mesh.parts(); p++) {
		std::stringstream ssDir;
		ssDir << ss.str() << "/" << p + _mesh.parts() * config::env::MPIrank;
		system(ssDir.str().c_str());
	}

	coordinates(_mesh.coordinates());
	elements(_mesh);
	materials(_mesh, _mesh.materials());
	boundaryConditions(_mesh.coordinates(), _mesh.boundaryConditions(), _mesh.DOFs());
	boundaries(_mesh);
}

void Esdata::coordinates(const Coordinates &coordinates)
{
	cilk_for (size_t p = 0; p < coordinates.parts(); p++) {
		std::ofstream os;
		eslocal value;
		esglobal index;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * config::env::MPIrank << "/coordinates.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		value = coordinates.localSize(p);
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		for (eslocal i = 0; i < coordinates.localSize(p); i++) {
			index = coordinates.globalIndex(i, p);
			os.write(reinterpret_cast<const char*>(&index), sizeof(esglobal));
			const Point &point = coordinates.get(i, p);
			os.write(reinterpret_cast<const char*>(&point), Point::size() * sizeof(double));
		}
		os.close();
	}
}


void Esdata::boundaryConditions(const Coordinates &coordinates, const std::vector<BoundaryCondition*> &conditions, size_t DOFs)
{
	cilk_for (size_t p = 0; p < coordinates.parts(); p++) {
		std::ofstream os;
		eslocal value;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * config::env::MPIrank << "/boundaryConditions.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		ESINFO(GLOBAL_ERROR) << "Broken VTK OUTPUT";
//		size_t counter = conditions.size();
//		os.write(reinterpret_cast<const char*>(&counter), sizeof(size_t));
//		auto &l2c = coordinates.localToCluster(p);
//		for (size_t i = 0; i < conditions.size(); i++) {
//			if (conditions[i]->type() == ConditionType::DIRICHLET) {
//				eslocal DOFs = 3;
//				size_t size = 0;
//				for (size_t j = 0; j < conditions[i]->DOFs().size(); j++) {
//					if (std::binary_search(l2c.begin(), l2c.end(), conditions[i]->DOFs()[j] / DOFs)) {
//						size++;
//					}
//				}
//
//				os.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
//				double dirichletValue = conditions[i]->value(Point3D());
//				os.write(reinterpret_cast<const char*>(&dirichletValue), sizeof(double));
//
//				for (size_t j = 0; j < conditions[i]->DOFs().size(); j++) {
//					if (std::binary_search(l2c.begin(), l2c.end(), conditions[i]->DOFs()[j] / DOFs)) {
//						value = std::lower_bound(l2c.begin(), l2c.end(), conditions[i]->DOFs()[j] / DOFs) - l2c.begin();
//						value = DOFs * value + conditions[i]->DOFs()[j] % DOFs;
//						os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
//					}
//				}
//			}
//		}

		os.close();
	}
}

void Esdata::elements(const Mesh &mesh)
{
	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		const std::vector<eslocal> &parts = mesh.getPartition();
		const std::vector<Element*> &elements = mesh.getElements();
		eslocal value;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * config::env::MPIrank << "/elements.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		value = parts[p + 1] - parts[p];
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			os << *(elements[e]);
		}
		os.close();
	}
}

void Esdata::materials(const Mesh &mesh, const std::vector<Material> &materials)
{
	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * config::env::MPIrank << "/materials.dat";
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		int size = materials.size();
		os.write(reinterpret_cast<const char*>(&size), sizeof(int));
		for (size_t i = 0; i < materials.size(); i++) {
			os.write(reinterpret_cast<const char*>(&materials[i].density), sizeof(double));
			os.write(reinterpret_cast<const char*>(&materials[i].youngModulus), sizeof(double));
			os.write(reinterpret_cast<const char*>(&materials[i].poissonRatio), sizeof(double));
			os.write(reinterpret_cast<const char*>(&materials[i].termalExpansion), sizeof(double));
			os.write(reinterpret_cast<const char*>(&materials[i].termalCapacity), sizeof(double));
			os.write(reinterpret_cast<const char*>(&materials[i].termalConduction), sizeof(double));
		}
		os.close();
	}
}

void Esdata::boundaries(const Mesh &mesh)
{
	Boundaries boundaries = mesh.subdomainBoundaries();
	const Boundaries &cBoundaries = mesh.clusterBoundaries();

	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, boundaries.size());

	std::vector<std::vector<std::vector<esglobal> > > sBuffer(threads, std::vector<std::vector<esglobal> > (mesh.neighbours().size()));
	std::vector<std::vector<esglobal> > rBuffer(mesh.neighbours().size());

	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			for (size_t j = 0; j < boundaries[i].size(); j++) {
				boundaries[i][j] += config::env::MPIrank * mesh.parts();
			}
			for (size_t n = 0; n < mesh.neighbours().size(); n++) {
				if (std::binary_search(cBoundaries[i].begin(), cBoundaries[i].end(), mesh.neighbours()[n])) {
					sBuffer[t][n].push_back(mesh.coordinates().globalIndex(i));
					sBuffer[t][n].push_back(boundaries[i].size());
					sBuffer[t][n].insert(sBuffer[t][n].end(), boundaries[i].begin(), boundaries[i].end());
				}
			}
		}
	}

	cilk_for (size_t n = 0; n < mesh.neighbours().size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	std::vector<MPI_Request> req(mesh.neighbours().size());

	for (size_t n = 0; n < mesh.neighbours().size(); n++) {
		MPI_Isend(sBuffer[0][n].data(), sizeof(esglobal) * sBuffer[0][n].size(), MPI_BYTE, mesh.neighbours()[n], 0, MPI_COMM_WORLD, req.data() + n);
	}

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(mesh.neighbours().begin(), mesh.neighbours().end(), neighbour) - mesh.neighbours().begin();
	};

	int flag, counter = 0;
	MPI_Status status;
	while (counter < mesh.neighbours().size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(esglobal));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(mesh.neighbours().size(), req.data(), MPI_STATUSES_IGNORE);

	for (size_t n = 0; n < mesh.neighbours().size(); n++) {
		for (size_t i = 0; i < rBuffer[n].size(); i++) {
			esglobal index = mesh.coordinates().clusterIndex(rBuffer[n][i]);
			boundaries[index].insert(boundaries[index].end(), rBuffer[n].begin() + i + 2, rBuffer[n].begin() + i + 2 + rBuffer[n][i + 1]);
			i += rBuffer[n][i + 1] + 1;
		}
	}

	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			std::sort(boundaries[i].begin(), boundaries[i].end());
		}
	}

	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		eslocal value, size;
		esglobal index;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * config::env::MPIrank << "/clusterBoundaries.dat";
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		size = 0;
		for (size_t i = 0; i < boundaries.size(); i++) {
			if (std::binary_search(boundaries[i].begin(), boundaries[i].end(), p + _mesh.parts() * config::env::MPIrank)) {
				size++;
			}
		}
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));

		for (size_t i = 0; i < boundaries.size(); i++) {
			if (std::binary_search(boundaries[i].begin(), boundaries[i].end(), p + _mesh.parts() * config::env::MPIrank)) {
				size = boundaries[i].size();
				os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
				for (auto it = boundaries[i].begin(); it != boundaries[i].end(); ++it) {
					value = *it;
					os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
				}
			}
		}
		os.close();
	}
}

