
#include "esdata.h"

using namespace esoutput;

void Esdata::store(const mesh::Mesh &mesh)
{
	std::stringstream ss;
	ss << "mkdir -p " << _path;

	system(ss.str().c_str());

	eslocal value;
	esglobal index;

	for (size_t p = 0; p < mesh.parts(); p++) {
		std::stringstream ssDir;
		ssDir << ss.str() << "/" << p;
		system(ssDir.str().c_str());

		std::stringstream ssMesh;
		ssMesh << _path << "/" << p << "/elements.dat";
		std::stringstream ssPoints;
		ssPoints << _path << "/" << p << "/coordinates.dat";
		std::ofstream os;

		// save elements
		const std::vector<eslocal> &parts = mesh.getPartition();
		const std::vector<mesh::Element*> &elements = mesh.getElements();
		value = parts[p + 1] - parts[p];

		os.open(ssMesh.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			os << *(elements[e]);
		}
		os.close();

		// save coordinates
		const std::vector<eslocal> &l2c = mesh.coordinates().localToCluster(p);
		value = mesh.coordinates().localSize(p);

		os.open(ssPoints.str().c_str(), std::ofstream::binary | std::ofstream::trunc);
		os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
		for (eslocal i = 0; i < l2c.size(); i++) {
			index = mesh.coordinates().globalIndex(i, p);
			os.write(reinterpret_cast<const char*>(&index), sizeof(esglobal));
			const mesh::Point &point = mesh.coordinates().get(i, p);
			os.write(reinterpret_cast<const char*>(&point), mesh::Point::size() * sizeof(double));
		}
		os.close();

//		// save coordinates' properties
//		for (size_t i = 0; i < _coordinates.propertiesSize(); i++) {
//			const std::map<eslocal, double> &property = _coordinates.property(i).values();
//			eslocal counter = 0;
//			const std::vector<eslocal> &l2c = _coordinates.localToCluster(p);
//			for (size_t j = 0; j < l2c.size(); j++) {
//				if (property.find(l2c[j]) != property.end()) {
//					counter++;
//				}
//			}
//			os.write(reinterpret_cast<const char*>(&counter), sizeof(eslocal));
//			for (eslocal j = 0; j < l2c.size(); j++) {
//				if (property.find(l2c[j]) != property.end()) {
//					os.write(reinterpret_cast<const char*>(&j), sizeof(eslocal));
//					os.write(reinterpret_cast<const char*>(&property.find(l2c[j])->second), sizeof(double));
//				}
//			}
//		}
	}
}

