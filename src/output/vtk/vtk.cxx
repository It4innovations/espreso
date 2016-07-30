
// Dummy VTK file
// Always store VTK Legacy format

#include "vtk.h"

using namespace espreso::output;

static espreso::Point clusterCenter(const espreso::Coordinates &coordinates)
{
	espreso::Point center(0, 0, 0);

	for (size_t i = 0; i < coordinates.clusterSize(); i++) {
		center += coordinates[i];
	}
	center /= coordinates.clusterSize();

	return center;
}

static std::vector<espreso::Point> subdomainsCenter(const espreso::Coordinates &coordinates)
{
	std::vector<espreso::Point> centers(coordinates.parts(), espreso::Point(0, 0, 0));

	for (size_t p = 0; p < coordinates.parts(); p++) {
		for (size_t i = 0; i < coordinates.localSize(p); i++) {
			centers[p] += coordinates.get(i, p);
		}
		centers[p] /= coordinates.localSize(p);
	}

	return centers;
}

static espreso::Point shrink(
		const espreso::Point &p,
		const espreso::Point &subdomainCenter, double subdomainShrinkRatio,
		const espreso::Point &clusterCenter, double clusterShrinkRatio)
{
	espreso::Point x = p;
	x = subdomainCenter + (x - subdomainCenter) * subdomainShrinkRatio;
	x = clusterCenter + (x - clusterCenter) * clusterShrinkRatio;
	return x;
}

static void head(std::ofstream &os)
{
	os << "# vtk DataFile Version 3.0\n";
	os << "Test\n";
	os << "ASCII\n";
	os << "\n";
}

static void coordinates(std::ofstream &os, const espreso::Coordinates &coordinates, double shrinkSubdomain, double shringCluster)
{
	size_t parts = coordinates.parts();

	size_t cSize = 0;
	for (size_t p = 0; p < parts; p++) {
		cSize += coordinates.localToCluster(p).size();
	}

	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << cSize << " float\n";

	espreso::Point cCenter = clusterCenter(coordinates);
	std::vector<espreso::Point> sCenters = subdomainsCenter(coordinates);

	for (size_t p = 0; p < parts; p++) {
		for (size_t i = 0; i < coordinates.localToCluster(p).size(); i++) {
			os << shrink(coordinates.get(i, p), sCenters[p], shrinkSubdomain, cCenter, shringCluster) << "\n";
		}
	}
	os << "\n";
}

static void coordinates(std::ofstream &os, const espreso::Coordinates &coordinates, const std::vector<std::vector<eslocal> > &points, double shrinkSubdomain, double shringCluster)
{
	size_t parts = coordinates.parts();

	size_t cSize = 0;
	for (size_t p = 0; p < parts; p++) {
		cSize += points[p].size();
	}

	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << cSize << " float\n";

	espreso::Point cCenter = clusterCenter(coordinates);
	std::vector<espreso::Point> sCenters = subdomainsCenter(coordinates);

	for (size_t p = 0; p < parts; p++) {
		for (size_t i = 0; i < points[p].size(); i++) {
			os << shrink(coordinates.get(points[p][i], p), sCenters[p], shrinkSubdomain, cCenter, shringCluster) << "\n";
		}
	}
	os << "\n";
}

static void points(std::ofstream &os, const std::vector<std::vector<eslocal> > &points)
{
	size_t size = 0;
	for (size_t p = 0; p < points.size(); p++) {
		size += points[p].size() + 1;
	}

	size_t offset = 0;
	os << "CELLS " << points.size() << " " << size << "\n";
	for (size_t p = 0; p < points.size(); p++) {
		os << points[p].size();
		for (size_t j = 0; j < points[p].size(); j++) {
			os << " " << offset + j;
		}
		os << "\n";

		offset += points[p].size();
	}

	os << "\n";

	os << "CELL_TYPES " << points.size() << "\n";
	for (size_t p = 0; p < points.size(); p++) {
		os << "2\n";
	}
	os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	os << "CELL_DATA " << points.size() << "\n";
	os << "SCALARS decomposition int 1\n";
	os << "LOOKUP_TABLE decomposition\n";
	for (size_t p = 0; p < points.size(); p++) {
		os << p << "\n";
	}
	os << "\n";
}

static void elements(std::ofstream &os, const espreso::Mesh &mesh)
{
	const std::vector<espreso::Element*> &elements = mesh.getElements();
	const std::vector<eslocal> &partition = mesh.getPartition();
	size_t parts = mesh.parts();

	size_t size = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		size += elements[i]->size() + 1;
	}

	// ELEMENTS
	size_t offset = 0;
	os << "CELLS " << elements.size() << " " << size << "\n";
	for (size_t p = 0; p < parts; p++) {
		for (size_t i = partition[p]; i < partition[p + 1]; i++) {
			// elements
			os << elements[i]->size();
			for (size_t j = 0; j < elements[i]->size(); j++) {
				os << " " << elements[i]->node(j) + offset;
			}
			os << "\n";
		}

		offset += mesh.coordinates().localSize(p);
	}

	os << "\n";

	// ELEMENTS TYPES
	os << "CELL_TYPES " << elements.size() << "\n";
	for (size_t p = 0; p < parts; p++) {
		for (size_t i = partition[p]; i < partition[p + 1]; i++) {
			os << elements[i]->vtkCode() << "\n";
		}
	}
	os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	os << "CELL_DATA " << elements.size() << "\n";
	os << "SCALARS decomposition int 1\n";
	os << "LOOKUP_TABLE decomposition\n";
	for (size_t p = 0; p < parts; p++) {
		for (eslocal i = partition[p]; i < partition[p + 1]; i++) {
			os << p << "\n";

		}
	}
	os << "\n";

	// DECOMPOSITION TO MATERIAL
	os << "SCALARS materials int 1\n";
	os << "LOOKUP_TABLE materials\n";
	for (size_t p = 0; p < parts; p++) {
		for (eslocal i = partition[p]; i < partition[p + 1]; i++) {
			os << elements[i]->getParam(espreso::Element::MATERIAL) << "\n";

		}
	}
	os << "\n";
}

static void coordinatesDisplacement(std::ofstream &os, const std::vector<std::vector<double> > &displacement, size_t DOFs)
{
	size_t size = 0;
	for (size_t p = 0; p < displacement.size(); p++) {
		size += displacement[p].size() / DOFs;
	}

	os << "\n";
	os << "POINT_DATA " << size << "\n";
	os << "SCALARS displacements float " << DOFs << "\n";
	os << "LOOKUP_TABLE default\n";
	for (size_t p = 0; p < displacement.size(); p++) {
		for (size_t i = 0; i < displacement[p].size() / DOFs; i++) {
			for (size_t d = 0; d < DOFs; d++) {
				os << displacement[p][DOFs * i + d] << " ";
			}
			os << "\n";
		}
	}
}

VTK::VTK(const Mesh &mesh, const std::string &path): Results(mesh, path)
{
	switch (config::output::OUTPUT_FORMAT) {
	case config::output::OUTPUT_FORMATAlternatives::VTK_LEGACY_FORMAT:
		break;
	default:
		ESINFO(ALWAYS) << TextColor::YELLOW << "Warning: ESPRESO not contains a library for saving generic VTK format. VTK Legacy format is used.";
	}
}

void VTK::store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster)
{
	std::stringstream ss;
	ss << _path << config::env::MPIrank;
	if (config::solver::TIME_STEPS > 1) {
		ss << "_" << counter++;
	}
	ss << ".vtk";

	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head(os);
	coordinates(os, _mesh.coordinates(),shrinkSubdomain, shrinkCluster);
	elements(os, _mesh);
	coordinatesDisplacement(os, displacement, displacement[0].size() / _mesh.coordinates().localSize(0));

	os.close();
}

void VTK::mesh(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shrinkCluster)
{
	std::stringstream ss;
	ss << path << config::env::MPIrank << ".vtk";

	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head(os);
	coordinates(os, mesh.coordinates(),shrinkSubdomain, shrinkCluster);
	elements(os, mesh);

	os.close();
}

void VTK::properties(const Mesh &mesh, const std::string &path, std::vector<Property> properties, double shrinkSubdomain, double shrinkCluster)
{
	std::stringstream ss;
	ss << path << config::env::MPIrank << ".vtk";

	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head(os);
	coordinates(os, mesh.coordinates(),shrinkSubdomain, shrinkCluster);
	elements(os, mesh);

	// TODO:

	os.close();
}

void VTK::fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	std::vector<std::vector<eslocal> > fixPoints(mesh.parts());
	for (size_t p = 0; p < mesh.parts(); p++) {
		fixPoints[p] = mesh.computeFixPoints(p, config::mesh::FIX_POINTS);
	}

	std::stringstream ss;
	ss << path << config::env::MPIrank << ".vtk";
	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head(os);
	coordinates(os, mesh.coordinates(), fixPoints, shrinkSubdomain, shringCluster);
	points(os, fixPoints);

	os.close();
}

void VTK::corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	std::vector<std::vector<eslocal> > corners(mesh.parts());
	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t i = 0; i < mesh.coordinates().localToCluster(p).size(); i++) {
			if (mesh.subdomainBoundaries().isCorner(mesh.coordinates().localToCluster(p)[i])) {
				corners[p].push_back(i);
			}
		}
	}

	std::stringstream ss;
	ss << path << config::env::MPIrank << ".vtk";
	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head(os);
	coordinates(os, mesh.coordinates(), corners, shrinkSubdomain, shringCluster);
	points(os, corners);

	os.close();
}
