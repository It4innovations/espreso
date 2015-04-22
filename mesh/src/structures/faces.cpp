#include "faces.h"

Faces::Faces(const Mesh &mesh, const Coordinates &coordinates): _coordinates(coordinates)
{
	const std::vector<idx_t> &parts = mesh.getPartition();
	const std::vector<Element*> &elements = mesh.getElements();
	for (size_t p = 0; p < parts.size() - 1; p++) {
		for (idx_t i = parts[p]; i < parts[p + 1]; i++) {
			elements[i]->fillFaces(_faces, p);
		}
	}
}

Faces::Faces(const Mesh &mesh, const Coordinates &coordinates, const BoundaryNodes &nodes): _coordinates(coordinates)
{
	const std::vector<idx_t> &parts = mesh.getPartition();
	const std::vector<Element*> &elements = mesh.getElements();
	for (size_t p = 0; p < parts.size() - 1; p++) {
		for (idx_t i = parts[p]; i < parts[p + 1]; i++) {
			elements[i]->fillFacesOnBorder(_faces, nodes, p);
		}
	}
	BoundaryFaces::iterator it = _faces.begin();
	while (it != _faces.end()) {
		if (it->second.first == it->second.second || it->second.second == - 1) {
			BoundaryFaces::iterator rm = it;
			++it;
			_faces.erase(rm);
		} else {
			++it;
		}
	}
}

Faces::~Faces()
{
	BoundaryFaces::iterator it;
	for (it = _faces.begin(); it != _faces.end(); ++it) {
		delete it->first;
	}
}

std::ostream& operator<<(std::ostream& os, const Faces &f)
{
	size_t inner = 0, border = 0, outer = 0;
	BoundaryFaces::const_iterator it;
	for (it = f._faces.begin(); it != f._faces.end(); ++it) {
		if (it->second.second == -1) {
			outer++;
		} else if (it->second.first != it->second.second) {
			border++;
		}
		if (it->second.first == it->second.second) {
			inner++;
		}
		os << *(it->first) << " " << it->second.first << " " << it->second.second << "\n";
	}
	os << "inner: " << inner << "\n";
	os << "border: " << border << "\n";
	os << "outer: " << outer << "\n";
	return os;
}



