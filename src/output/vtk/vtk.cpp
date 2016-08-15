
#include "vtk.h"

using namespace espreso::output;

VTK::VTK(const Mesh &mesh, const std::string &path): ResultStore(mesh, path)
{
	const Coordinates &coordinates = mesh.coordinates();

	for (size_t i = 0; i < coordinates.clusterSize(); i++) {
		_clusterCenter += coordinates[i];
	}
	_clusterCenter /= coordinates.clusterSize();

	_subdomainsCenter.resize(_mesh.parts());
	for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t i = 0; i < coordinates.localSize(p); i++) {
			_subdomainsCenter[p] += coordinates.get(i, p);
		}
		_subdomainsCenter[p] /= coordinates.localSize(p);
	}
}

espreso::Point VTK::shrink(const Point &p,
			const Point &subdomainCenter, double subdomainShrinkRatio,
			const Point &clusterCenter, double clusterShrinkRatio)
{
	Point x = p;
	x = subdomainCenter + (x - subdomainCenter) * subdomainShrinkRatio;
	x = clusterCenter + (x - clusterCenter) * clusterShrinkRatio;
	return x;
}

void VTK::store(const std::vector<std::vector<eslocal> > &points, double shrinkSubdomain, double shringCluster)
{
	std::stringstream ss;
	ss << _path << config::env::MPIrank << ".vtk";
	_vtk.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head();
	coordinates(_mesh.coordinates(), points, shrinkSubdomain, shringCluster);
	this->points(points);

	_vtk.close();
}


void VTK::store(double shrinkSubdomain, double shringCluster)
{
	std::stringstream ss;
	ss << _path << config::env::MPIrank << ".vtk";
	_vtk.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head();
	coordinates(_mesh.coordinates(), shrinkSubdomain, shringCluster);
	elements(_mesh);

	_vtk.close();
}


void VTK::store(std::vector<std::vector<double> > &displacement, size_t dofs, double shrinkSubdomain, double shringCluster)
{
	std::stringstream ss;
	ss << _path << config::env::MPIrank << ".vtk";
	_vtk.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head();
	coordinates(_mesh.coordinates(), shrinkSubdomain, shringCluster);
	elements(_mesh);
	coordinatesDisplacement(displacement, dofs);

	_vtk.close();
}

void VTK::head()
{
	_vtk << "# vtk DataFile Version 3.0\n";
	_vtk << "Test\n";
	_vtk << "ASCII\n";
	_vtk << "\n";
}

void VTK::coordinates(const Coordinates &coordinates, double shrinkSubdomain, double shringCluster)
{
	size_t parts = coordinates.parts();

	size_t cSize = 0;
	for (size_t p = 0; p < parts; p++) {
		cSize += coordinates.localToCluster(p).size();
	}

	_vtk << "DATASET UNSTRUCTURED_GRID\n";
	_vtk << "POINTS " << cSize << " float\n";

	for (size_t p = 0; p < parts; p++) {
		for (size_t i = 0; i < coordinates.localToCluster(p).size(); i++) {
			_vtk << shrink(coordinates.get(i, p), _subdomainsCenter[p], shrinkSubdomain, _clusterCenter, shringCluster) << "\n";
		}
	}
	_vtk << "\n";
}

void VTK::coordinates(const Coordinates &coordinates, const std::vector<std::vector<eslocal> > &points, double shrinkSubdomain, double shringCluster)
{
	size_t parts = coordinates.parts();

	size_t cSize = 0;
	for (size_t p = 0; p < parts; p++) {
		cSize += points[p].size();
	}

	_vtk << "DATASET UNSTRUCTURED_GRID\n";
	_vtk << "POINTS " << cSize << " float\n";

	for (size_t p = 0; p < parts; p++) {
		for (size_t i = 0; i < points[p].size(); i++) {
			_vtk << shrink(coordinates.get(points[p][i], p), _subdomainsCenter[p], shrinkSubdomain, _clusterCenter, shringCluster) << "\n";
		}
	}
	_vtk << "\n";
}

void VTK::elements(const Mesh &mesh)
{
	const std::vector<Element*> &elements = mesh.elements();
	const std::vector<eslocal> &partition = mesh.getPartition();
	size_t parts = mesh.parts();

	size_t size = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		size += elements[i]->nodes() + 1;
	}

	// ELEMENTS
	size_t offset = 0;
	_vtk << "CELLS " << elements.size() << " " << size << "\n";
	for (size_t p = 0; p < parts; p++) {
		for (size_t i = partition[p]; i < partition[p + 1]; i++) {
			// elements
			_vtk << elements[i]->nodes();
			for (size_t j = 0; j < elements[i]->nodes(); j++) {
				_vtk << " " << mesh.coordinates().localIndex(elements[i]->node(j), p) + offset;
			}
			_vtk << "\n";
		}

		offset += mesh.coordinates().localSize(p);
	}

	_vtk << "\n";

	// ELEMENTS TYPES
	_vtk << "CELL_TYPES " << elements.size() << "\n";
	for (size_t p = 0; p < parts; p++) {
		for (size_t i = partition[p]; i < partition[p + 1]; i++) {
			_vtk << elements[i]->vtkCode() << "\n";
		}
	}
	_vtk << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	_vtk << "CELL_DATA " << elements.size() << "\n";
	_vtk << "SCALARS decomposition int 1\n";
	_vtk << "LOOKUP_TABLE decomposition\n";
	for (size_t p = 0; p < parts; p++) {
		for (eslocal i = partition[p]; i < partition[p + 1]; i++) {
			_vtk << p << "\n";

		}
	}
	_vtk << "\n";

	// DECOMPOSITION TO MATERIAL
	_vtk << "SCALARS materials int 1\n";
	_vtk << "LOOKUP_TABLE materials\n";
	for (size_t p = 0; p < parts; p++) {
		for (eslocal i = partition[p]; i < partition[p + 1]; i++) {
			_vtk << elements[i]->param(Element::MATERIAL) << "\n";

		}
	}
	_vtk << "\n";
}

void VTK::points(const std::vector<std::vector<eslocal> > &points)
{
	size_t size = 0;
	for (size_t p = 0; p < points.size(); p++) {
		size += points[p].size() + 1;
	}

	size_t offset = 0;
	_vtk << "CELLS " << points.size() << " " << size << "\n";
	for (size_t p = 0; p < points.size(); p++) {
		_vtk << points[p].size();
		for (size_t j = 0; j < points[p].size(); j++) {
			_vtk << " " << offset + j;
		}
		_vtk << "\n";

		offset += points[p].size();
	}

	_vtk << "\n";

	_vtk << "CELL_TYPES " << points.size() << "\n";
	for (size_t p = 0; p < points.size(); p++) {
		_vtk << "2\n";
	}
	_vtk << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	_vtk << "CELL_DATA " << points.size() << "\n";
	_vtk << "SCALARS decomposition int 1\n";
	_vtk << "LOOKUP_TABLE decomposition\n";
	for (size_t p = 0; p < points.size(); p++) {
			_vtk << p << "\n";
	}
	_vtk << "\n";
}




