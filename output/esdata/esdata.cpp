
#include "esdata.h"

using namespace esoutput;

void Esdata::store(const mesh::Mesh &mesh, double shrinkSubdomain, double shrinkCluster)
{
	// shrink is ignored!!!

	std::stringstream ss;
	ss << "mkdir -p " << _path;
	system(ss.str().c_str());
	for (size_t p = 0; p < mesh.parts(); p++) {
		std::stringstream ssDir;
		ssDir << ss.str() << "/" << p;
		system(ssDir.str().c_str());
	}

	coordinates(mesh.coordinates());
	elements(mesh);
	boundaryConditions(mesh.coordinates());
	boundaries(mesh);
}

void Esdata::coordinates(const mesh::Coordinates &coordinates)
{
	std::ofstream os;
	eslocal value;
	esglobal index;

	for (size_t p = 0; p < coordinates.parts(); p++) {
		std::stringstream ss;
		ss << _path << "/" << p << "/coordinates.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		value = coordinates.localSize(p);
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		for (eslocal i = 0; i < coordinates.localSize(p); i++) {
			index = coordinates.globalIndex(i, p);
			os.write(reinterpret_cast<const char*>(&index), sizeof(esglobal));
			const mesh::Point &point = coordinates.get(i, p);
			os.write(reinterpret_cast<const char*>(&point), mesh::Point::size() * sizeof(double));
		}
		os.close();
	}
}


void Esdata::boundaryConditions(const mesh::Coordinates &coordinates)
{
	std::ofstream os;
	eslocal value;

	for (size_t p = 0; p < coordinates.parts(); p++) {
		std::stringstream ss;
		ss << _path << "/" << p << "/boundaryConditions.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		for (size_t i = 0; i < coordinates.propertiesSize(); i++) {
			const std::map<eslocal, double> &property = coordinates.property(i).values();
			eslocal counter = 0;
			for (size_t j = 0; j < coordinates.localSize(p); j++) {
				if (property.find(coordinates.clusterIndex(j, p)) != property.end()) {
					counter++;
				}
			}
			os.write(reinterpret_cast<const char*>(&counter), sizeof(eslocal));
			for (eslocal j = 0; j < coordinates.localSize(p); j++) {
				if (property.find(coordinates.clusterIndex(j, p)) != property.end()) {
					os.write(reinterpret_cast<const char*>(&j), sizeof(eslocal));
					os.write(reinterpret_cast<const char*>(&property.find(coordinates.clusterIndex(j, p))->second), sizeof(double));
				}
			}
		}
		os.close();
	}
}

void Esdata::elements(const mesh::Mesh &mesh)
{
	std::ofstream os;
	const std::vector<eslocal> &parts = mesh.getPartition();
	const std::vector<mesh::Element*> &elements = mesh.getElements();
	eslocal value;

	for (size_t p = 0; p < mesh.parts(); p++) {
		std::stringstream ss;
		ss << _path << "/" << p << "/elements.dat";

		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		value = parts[p + 1] - parts[p];
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			os << *(elements[e]);
		}
		os.close();
	}
}

void Esdata::boundaries(const mesh::Mesh &mesh)
{
	std::ofstream os;
	const mesh::Boundaries &boundaries = mesh.subdomainBoundaries();
	eslocal value, size;
	esglobal index;

	for (size_t p = 0; p < mesh.parts(); p++) {
		std::stringstream ss;
		ss << _path << "/" << p << "/clusterBoundaries.dat";
		std::cout << "X " << p << " " << ss.str() << std::endl;
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		size = 0;
		for (size_t i = 0; i < boundaries.size(); i++) {
			if (boundaries[i].find(p) != boundaries[i].end()) {
				size++;
			}
		}
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));

		std::set<eslocal>::const_iterator it;
		for (size_t i = 0; i < boundaries.size(); i++) {
			if (boundaries[i].find(p) != boundaries[i].end()) {
				size = boundaries[i].size();
				os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
				for (it = boundaries[i].begin(); it != boundaries[i].end(); ++it) {
					value = *it;
					os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
				}
			}
		}
		os.close();
	}
}

