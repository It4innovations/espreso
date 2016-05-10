
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
		ssDir << ss.str() << "/" << p;
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
		ss << _path << "/" << p << "/coordinates.dat";

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

		size_t counter = conditions.size();
		os.write(reinterpret_cast<const char*>(&counter), sizeof(size_t));
		auto &l2c = coordinates.localToCluster(p);
		for (size_t i = 0; i < conditions.size(); i++) {
			if (conditions[i]->type() == ConditionType::DIRICHLET) {
				size_t size = 0;
				for (size_t j = 0; j < conditions[i]->DOFs().size(); j++) {
					if (std::binary_search(l2c.begin(), l2c.end(), conditions[i]->DOFs()[j] / DOFs)) {
						size++;
					}
				}

				os.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
				double dirichletValue = conditions[i]->value();
				os.write(reinterpret_cast<const char*>(&dirichletValue), sizeof(double));

				for (size_t j = 0; j < conditions[i]->DOFs().size(); j++) {
					if (std::binary_search(l2c.begin(), l2c.end(), conditions[i]->DOFs()[j] / DOFs)) {
						value = std::lower_bound(l2c.begin(), l2c.end(), conditions[i]->DOFs()[j] / DOFs) - l2c.begin();
						value = DOFs * value + conditions[i]->DOFs()[j] % DOFs;
						os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
					}
				}
			}
		}

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

void Esdata::materials(const Mesh &mesh, const std::vector<Material> &materials)
{
	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		std::stringstream ss;
		ss << _path << "/" << p << "/materials.dat";
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
	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
		std::ofstream os;
		const Boundaries &boundaries = mesh.subdomainBoundaries();
		eslocal value, size;
		esglobal index;

		std::stringstream ss;
		ss << _path << "/" << p << "/clusterBoundaries.dat";
		os.open(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		size = 0;
		for (size_t i = 0; i < boundaries.size(); i++) {
			if (std::binary_search(boundaries[i].begin(), boundaries[i].end(), p)) {
				size++;
			}
		}
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));

		for (size_t i = 0; i < boundaries.size(); i++) {
			if (std::binary_search(boundaries[i].begin(), boundaries[i].end(), p)) {
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

