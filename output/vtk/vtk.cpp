
#include "vtk.h"

using namespace esoutput;

void VTK::store(const mesh::Mesh &mesh, double shrinkSubdomain, double shringCluster)
{
	std::stringstream ss;
	ss << _file << _rank << ".vtk";
	_vtk.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head();
	coordinates(mesh.coordinates(), shrinkSubdomain, shringCluster);
	elements(mesh);

	_vtk.close();
}


void VTK::store(const mesh::Mesh &mesh, std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shringCluster)
{
	std::stringstream ss;
	ss << _file << _rank << ".vtk";
	_vtk.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head();
	coordinates(mesh.coordinates(), shrinkSubdomain, shringCluster);
	elements(mesh);
	coordinatesDisplacement(displacement);

	_vtk.close();
}

void VTK::head()
{
	_vtk << "# vtk DataFile Version 3.0\n";
	_vtk << "Test\n";
	_vtk << "ASCII\n";
	_vtk << "\n";
}

void VTK::coordinates(const mesh::Coordinates &coordinates, double shrinkSubdomain, double shringCluster)
{
	size_t parts = coordinates.parts();

	size_t cSize = 0;
	for (size_t p = 0; p < parts; p++) {
		cSize += coordinates.localToCluster(p).size();
	}

	_vtk << "DATASET UNSTRUCTURED_GRID\n";
	_vtk << "POINTS " << cSize << " float\n";

	mesh::Point cCenter;
	for (size_t i = 0; i < coordinates.size(); i++) {
		cCenter += coordinates[i];
	}
	cCenter /= coordinates.size();


	for (size_t p = 0; p < parts; p++) {
		mesh::Point sCenter;
		for (size_t i = 0; i < coordinates.localSize(p); i++) {
			sCenter += coordinates.get(i, p);
		}
		sCenter /= coordinates.localSize(p);

		for (size_t i = 0; i < coordinates.localSize(p); i++) {
			mesh::Point x = coordinates.get(i, p);
			x = sCenter + (x - sCenter) * shrinkSubdomain;
			x = cCenter + (x - cCenter) * shringCluster;
			_vtk << x << "\n";
		}
	}
	_vtk << "\n";
}

void VTK::elements(const mesh::Mesh &mesh)
{
	const std::vector<mesh::Element*> &elements = mesh.getElements();
	const mesh::Coordinates &coordinates = mesh.coordinates();
	const std::vector<eslocal> &partition = mesh.getPartition();
	const std::vector<eslocal> &fixPoints = mesh.getFixPoints();
	const mesh::Boundaries &boundaries = mesh.subdomainBoundaries();
	size_t parts = mesh.parts();

	size_t size = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		size += elements[i]->size() + 1;
	}
	size += fixPoints.size() + parts;

	size_t corners = 0;
	for (size_t i = 0; i < boundaries.size(); i++) {
		if (boundaries.isCorner(i)) {
			corners += boundaries[i].size();
		}
	}
	size += corners + 1;

	// ELEMENTS
	size_t offset = 0;
	_vtk << "CELLS " << elements.size() + parts + 1 << " " << size << "\n";
	for (size_t p = 0; p < parts; p++) {
		for (size_t i = partition[p]; i < partition[p + 1]; i++) {
			// elements
			_vtk << elements[i]->size();
			for (size_t j = 0; j < elements[i]->size(); j++) {
				_vtk << " " << elements[i]->node(j) + offset;
			}
			_vtk << "\n";
		}
		// fix points for part 'p'
		_vtk << fixPoints.size() / parts;
		for (size_t j = 0; j < fixPoints.size() / parts; j++) {
			_vtk << " " << fixPoints[p * fixPoints.size() / parts + j] + offset;
		}
		_vtk << "\n";

		offset += mesh.coordinates().localSize(p);
	}
	_vtk << corners;
	offset = 0;
	for (size_t p = 0; p < parts; p++) {
		for (size_t i = 0; i < boundaries.size(); i++) {
			if (boundaries.isCorner(i) && boundaries[i].count(p)) {
				_vtk << " " << coordinates.localIndex(i, p) + offset;
			}
		}
		offset += mesh.coordinates().localSize(p);
	}
	_vtk << "\n";

	_vtk << "\n";

	// ELEMENTS TYPES
	_vtk << "CELL_TYPES " << elements.size() + parts + 1 << "\n";
	for (size_t p = 0; p < parts; p++) {
		for (size_t i = partition[p]; i < partition[p + 1]; i++) {
			_vtk << elements[i]->vtkCode() << "\n";
		}
		_vtk << "2\n";
	}
	_vtk << "2\n";
	_vtk << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	_vtk << "CELL_DATA " << elements.size() + parts + 1 << "\n";
	_vtk << "SCALARS decomposition int 1\n";
	_vtk << "LOOKUP_TABLE decomposition\n";
	for (size_t p = 0; p < parts; p++) {
		for (eslocal i = partition[p]; i < partition[p + 1]; i++) {
			_vtk << p << "\n";

		}
		_vtk << parts << "\n";
	}
	_vtk << parts + 1 << "\n";
	_vtk << "\n";
}

void VTK::coordinatesDisplacement(const std::vector<std::vector<double> > &displacement)
{
	size_t size = 0;
	for (size_t p = 0; p < displacement.size(); p++) {
		size += displacement[p].size() / 3;
	}

	_vtk << "\n";
	_vtk << "POINT_DATA " << size << "\n";
	_vtk << "SCALARS displacements float 3\n";
	_vtk << "LOOKUP_TABLE default\n";
	for (size_t p = 0; p < displacement.size(); p++) {
		for (size_t i = 0; i < displacement[p].size() / 3; i++) {
			_vtk << displacement[p][3 * i + 0] << " ";
			_vtk << displacement[p][3 * i + 1] << " ";
			_vtk << displacement[p][3 * i + 2] << "\n";
		}
	}

//	size_t size = 0;
//	for (size_t p = 0; p < displacement.size(); p++) {
//		size += displacement[p].size() / 1;
//	}
//
//	_vtk << "\n";
//	_vtk << "POINT_DATA " << size << "\n";
//	_vtk << "SCALARS displacements float 1\n";
//	_vtk << "LOOKUP_TABLE default\n";
//	for (size_t p = 0; p < displacement.size(); p++) {
//		for (size_t i = 0; i < displacement[p].size(); i++) {
//			_vtk << displacement[p][i] << "\n";
//		}
//	}

	_vtk << "\n";
}




