#include "boundaries.h"

using namespace mesh;

Boundaries::Boundaries(Mesh &m)
	: _mesh(m), _boundaries(m.coordinates().size()), _corners(m.coordinates().size(), false)
{
	const std::vector<eslocal> &parts = m.getPartition();
	const std::vector<Element*> &elements = m.getElements();
	const Coordinates &c = m.coordinates();

	for (size_t p = 0; p + 1 < parts.size(); p++) {
		for (eslocal e = parts[p]; e < parts[p + 1]; e++) {
			for (size_t n = 0; n < elements[e]->size(); n++) {
				_boundaries[c.clusterIndex(elements[e]->node(n), p)].insert(p);
			}
		}
	}
}

void Boundaries::saveData()
{
	eslocal value, size;
	esglobal index;

	for (size_t p = 0; p < _mesh.parts(); p++) {
		std::stringstream ss;
		ss << "boundaries_" << p << ".dat";
		std::ofstream os(ss.str().c_str(), std::ofstream::binary | std::ofstream::trunc);

		size = 0;
		for (size_t i = 0; i < _boundaries.size(); i++) {
			if (_boundaries[i].find(p) != _boundaries[i].end()) {
				size++;
			}
		}
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));

		std::set<eslocal>::const_iterator it;
		for (size_t i = 0; i < _boundaries.size(); i++) {
			if (_boundaries[i].find(p) != _boundaries[i].end()) {
				size = _boundaries[i].size();
				os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
				for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
					value = *it;
					os.write(reinterpret_cast<const char*>(&value), sizeof(eslocal));
				}
			}
		}
	}
}

void Boundaries::loadData(const char *filename)
{
	_boundaries.clear();
	_corners.clear();

	std::ifstream is(filename, std::ifstream::binary);

	eslocal size, value;

	is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));

	_boundaries.resize(size);
	for (size_t i = 0; i < _boundaries.size(); i++) {
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal j = 0; j < size; j++) {
			is.read(reinterpret_cast<char *>(&value), sizeof(eslocal));
			_boundaries[i].insert(value);
		}
	}
}

std::ostream& mesh::operator<<(std::ostream& os, const Boundaries &b)
{
	std::set<eslocal>::const_iterator it;

	for (size_t i = 0; i < b._boundaries.size(); i++) {
		os << i << ": ";
		for (it = b._boundaries[i].begin(); it != b._boundaries[i].end(); ++it) {
			os << *it << " ";
		}
		os << "\n";
	}
	return os;
}




