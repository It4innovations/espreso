
#include "espresobinaryformat.h"

#include <fstream>

#include "../../configuration/environment.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/material.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/structures/elementtypes.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/settings/evaluator.h"

using namespace espreso::output;

void ESPRESOBinaryFormat::prepareDirectories(const std::string &path, size_t parts)
{
	std::stringstream ss;
	ss << "mkdir -p " << path;
	int out = system(ss.str().c_str());
	for (size_t p = 0; p < parts; p++) {
		std::stringstream ssDir;
		ssDir << ss.str() << "/" << p + parts * environment->MPIrank;
		out = system(ssDir.str().c_str());
		if (out) {
			ESINFO(ERROR) << "Cannot create output directory '" << path << "/" << p + parts * environment->MPIrank << "'";
		}
	}
}

void ESPRESOBinaryFormat::store(const Mesh &mesh, const std::string &path)
{
	ESPRESOBinaryFormat(mesh, path);
}

ESPRESOBinaryFormat::ESPRESOBinaryFormat(const Mesh &mesh, const std::string &path)
: _mesh(mesh), _path(path)
{
	metafile();
	coordinates();
	elements();
	materials();
	regions();
	boundaries();
}

void ESPRESOBinaryFormat::metafile()
{
	std::ofstream os;
	os.open(_path + "/description.txt", std::ofstream::trunc);
	os << _mesh.parts() * environment->MPIsize << "\n";
	os.close();
}

void ESPRESOBinaryFormat::coordinates()
{
	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		std::ofstream os;
		eslocal size;
		esglobal index;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/coordinates.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
		if (!os.is_open()) {
			ESINFO(ERROR) << "Cannot open file '" << ss.str() << "'";
		}

		size = _mesh.coordinates().localSize(p);
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (size_t i = 0; i < _mesh.coordinates().localSize(p); i++) {
			index = _mesh.coordinates().globalIndex(i, p);
			os.write(reinterpret_cast<const char*>(&index), sizeof(esglobal));
			const Point &point = _mesh.coordinates().get(i, p);
			os.write(reinterpret_cast<const char*>(&point), Point::size() * sizeof(double));
		}
		os.close();
	}
}

void ESPRESOBinaryFormat::elements()
{
	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		std::ofstream os;
		const std::vector<eslocal> &parts = _mesh.getPartition();
		const std::vector<Element*> &elements = _mesh.elements();
		eslocal size;

		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/elements.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
		if (!os.is_open()) {
			ESINFO(ERROR) << "Cannot open file '" << ss.str() << "'";
		}

		// elements
		size = parts[p + 1] - parts[p];
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			elements[e]->store(os, _mesh.coordinates(), p);
		}

		os.close();
	}
}

void ESPRESOBinaryFormat::materials()
{
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh.parts(); p++) {
		std::ofstream os;
		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/materials.dat";
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
		if (!os.is_open()) {
			ESINFO(ERROR) << "Cannot open file '" << ss.str() << "'";
		}

		eslocal size = _mesh.materials().size();
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (size_t i = 0; i < _mesh.materials().size(); i++) {
			_mesh.materials()[i]->store(os);
		}
		os.close();
	}
}

void ESPRESOBinaryFormat::regions()
{
	auto computeIntervals = [] (const std::vector<espreso::Element*> &elements, eslocal p) {
		std::vector<size_t> intervals;
		if (!elements.size()) {
			return intervals;
		}

		if (std::find(elements[0]->domains().begin(), elements[0]->domains().end(), p) != elements[0]->domains().end()) {
			intervals.push_back(0);
		}
		for (size_t i = 1; i < elements.size(); i++) {
			if (elements[i - 1]->domains() != elements[i]->domains()) {
				if (std::find(elements[i]->domains().begin(), elements[i]->domains().end(), p) != elements[i]->domains().end()) {
					if (intervals.size() % 2 == 0) {
						intervals.push_back(i);
					}
				} else {
					if (intervals.size() % 2 == 1) {
						intervals.push_back(i);
					}
				}
			}
		}
		if (intervals.size() % 2 == 1) {
			intervals.push_back(elements.size());
		}
		return intervals;
	};

	auto intervalsSize = [] (const std::vector<size_t> &intervals) {
		size_t size = 0;
		for (size_t i = 0; i < intervals.size(); i += 2) {
			size += intervals[2 * i + 1] - intervals[2 * i];
		}
		return size;
	};

	auto region2index = [&] (const Region *region) {
		auto it = std::find(_mesh.regions().begin(), _mesh.regions().end(), region);
		return it - _mesh.regions().begin();
	};

	auto evaluator2index = [&] (const Evaluator *evaluator) {
		auto it = std::find(_mesh.evaluators().begin(), _mesh.evaluators().end(), evaluator);
		return it - _mesh.evaluators().begin();
	};


	#pragma omp parallel for
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		std::ofstream os;
		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/regions.dat";
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
		if (!os.is_open()) {
			ESINFO(ERROR) << "Cannot open file '" << ss.str() << "'";
		}

		eslocal size;
		std::vector<Element*> faces;
		std::vector<Element*> edges;
		for (size_t i = 0; i < _mesh.faces().size(); i++) {
			if (_mesh.faces()[i]->regions().size()) {
				faces.push_back(_mesh.faces()[i]);
			}
		}
		for (size_t i = 0; i < _mesh.edges().size(); i++) {
			if (_mesh.edges()[i]->regions().size()) {
				edges.push_back(_mesh.edges()[i]);
			}
		}

		std::sort(faces.begin(), faces.end(), [] (Element *e1, Element *e2) { return e1->domains() < e2->domains(); });
		std::sort(edges.begin(), edges.end(), [] (Element *e1, Element *e2) { return e1->domains() < e2->domains(); });

		// faces
		std::vector<size_t> fIntervals = computeIntervals(faces, p);
		eslocal fSize = intervalsSize(fIntervals);
		os.write(reinterpret_cast<const char*>(&fSize), sizeof(eslocal));
		for (size_t i = 0; i < fIntervals.size(); i += 2) {
			for (size_t f = fIntervals[2 * i]; f < fIntervals[2 * i + 1]; f++) {
				faces[f]->store(os, _mesh.coordinates(), p);
			}
		}

		// edges
		std::vector<size_t> eIntervals = computeIntervals(edges, p);
		eslocal eSize = intervalsSize(eIntervals);
		os.write(reinterpret_cast<const char*>(&eSize), sizeof(eslocal));
		for (size_t i = 0; i < eIntervals.size(); i += 2) {
			for (size_t e = eIntervals[2 * i]; e < eIntervals[2 * i + 1]; e++) {
				edges[e]->store(os, _mesh.coordinates(), p);
			}
		}

		size = _mesh.evaluators().size();
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (size_t i = 0; i < _mesh.evaluators().size(); i++) {
			_mesh.evaluators()[i]->store(os);
		}

		// regions
		size = _mesh.regions().size();
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (size_t i = 0; i < _mesh.regions().size(); i++) {
			if (i > 1) {
				// not store default regions ELL_ELEMENT and ALL_NODES name
				eslocal length = _mesh.regions()[i]->name.size();
				os.write(reinterpret_cast<const char *>(&length), sizeof(eslocal));
				os.write(_mesh.regions()[i]->name.c_str(), _mesh.regions()[i]->name.size());
				os.write(reinterpret_cast<const char *>(&_mesh.regions()[i]->eType), sizeof(ElementType));
			}
			size = _mesh.regions()[i]->settings.size();
			os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
			for (size_t step = 0; step < _mesh.regions()[i]->settings.size(); step++) {
				size = _mesh.regions()[i]->settings[step].size();
				os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
				for (auto it = _mesh.regions()[i]->settings[step].begin(); it != _mesh.regions()[i]->settings[step].end(); ++it) {
					int property = (int)it->first;
					os.write(reinterpret_cast<const char *>(&property), sizeof(int));
					size = it->second.size();
					os.write(reinterpret_cast<const char *>(&size), sizeof(eslocal));
					for (size_t j = 0; j < it->second.size(); j++) {
						int index = evaluator2index(it->second[j]);
						os.write(reinterpret_cast<const char *>(&index), sizeof(int));
					}
				}
			}
		}

		// elements
		const std::vector<eslocal> &parts = _mesh.getPartition();
		const std::vector<Element*> &elements = _mesh.elements();
		size = parts[p + 1] - parts[p];
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			size = elements[e]->regions().size() - 1;
			os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
			for (size_t r = 1; r < elements[e]->regions().size(); r++) {
				int index = region2index(elements[e]->regions()[r]);
				os.write(reinterpret_cast<const char*>(&index), sizeof(int));
			}
		}

		// faces
		os.write(reinterpret_cast<const char*>(&fSize), sizeof(eslocal));
		for (size_t i = 0; i < fIntervals.size(); i += 2) {
			for (size_t f = fIntervals[2 * i]; f < fIntervals[2 * i + 1]; f++) {
				size = faces[f]->regions().size();
				os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
				for (size_t r = 0; r < faces[f]->regions().size(); r++) {
					int index = region2index(faces[f]->regions()[r]);
					os.write(reinterpret_cast<const char*>(&index), sizeof(int));
				}
			}
		}

		// edges
		os.write(reinterpret_cast<const char*>(&eSize), sizeof(eslocal));
		for (size_t i = 0; i < eIntervals.size(); i += 2) {
			for (size_t e = eIntervals[2 * i]; e < eIntervals[2 * i + 1]; e++) {
				size = edges[e]->regions().size();
				os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
				for (size_t r = 0; r < edges[e]->regions().size(); r++) {
					int index = region2index(edges[e]->regions()[r]);
					os.write(reinterpret_cast<const char*>(&index), sizeof(int));
				}
			}
		}

		// nodes
		size = _mesh.coordinates().localSize(p);
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (size_t i = 0; i < _mesh.coordinates().localSize(p); i++) {
			const Element* node = _mesh.nodes()[_mesh.coordinates().clusterIndex(i, p)];
			size = node->regions().size() - 1;
			os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
			for (size_t r = 1; r < node->regions().size(); r++) {
				int index = region2index(node->regions()[r]);
				os.write(reinterpret_cast<const char*>(&index), sizeof(int));
			}
		}

		os.close();
	}
}

void ESPRESOBinaryFormat::boundaries()
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
	for  (size_t p = 0; p < _mesh.parts(); p++) {
		std::ofstream os;
		std::stringstream ss;
		ss << _path << "/" << p + _mesh.parts() * environment->MPIrank << "/boundaries.dat";
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
		if (!os.is_open()) {
			ESINFO(ERROR) << "Cannot open file '" << ss.str() << "'";
		}

		eslocal domain, size;

		size = _mesh.coordinates().localSize(p);
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));

		for (size_t i = 0; i < _mesh.coordinates().localSize(p); i++) {
			const std::vector<eslocal> &domains = _mesh.nodes()[_mesh.coordinates().clusterIndex(i, p)]->domains();

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

