
// Dummy VTK file
// Always store VTK Legacy format

#include "../vtk.h"

using namespace espreso::output;

static void open(std::ofstream &os, const std::string &path, size_t timeStep)
{
	std::stringstream ss;
	ss << path << espreso::config::env::MPIrank;
	if (timeStep < -1) {
		ss << "_" << timeStep;
	}
	ss << ".vtk";

	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);
}

static void head(std::ofstream &os)
{
	os << "# vtk DataFile Version 4.0\n";
	os << "ESPRESO output\n";
	os << "ASCII\n";
	os << "\n";
}

static void coordinates(std::ofstream &os, const espreso::Coordinates &coordinates, std::function<espreso::Point(const espreso::Point&, size_t)> shrink)
{
	size_t parts = coordinates.parts();

	size_t cSize = 0;
	for (size_t p = 0; p < parts; p++) {
		cSize += coordinates.localToCluster(p).size();
	}

	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << cSize << " float\n";

	for (size_t p = 0; p < parts; p++) {
		for (size_t i = 0; i < coordinates.localToCluster(p).size(); i++) {
			os << shrink(coordinates.get(i, p), p) << "\n";
		}
	}
	os << "\n";
}

static void coordinates(std::ofstream &os, const espreso::Coordinates &coordinates, const std::vector<espreso::Element*> &nodes, std::function<espreso::Point(const espreso::Point&, size_t)> shrink)
{
	size_t parts = coordinates.parts();

	size_t cSize = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		cSize += nodes[i]->domains().size();
	}

	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << cSize << " float\n";

	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t d = 0; d < nodes[i]->domains().size(); d++) {
			os << shrink(coordinates[nodes[i]->node(0)], nodes[i]->domains()[d]) << "\n";
		}
	}
	os << "\n";
}

static void nodes(std::ofstream &os, const std::vector<espreso::Element*> &nodes, size_t domains)
{
	std::vector<size_t> domainSize(domains, 0);
	size_t size = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t d = 0; d < nodes[i]->domains().size(); d++) {
			domainSize[nodes[i]->domains()[d]] += nodes[i]->domains().size();
		}
	}

	for (size_t i = 0; i < domains; i++) {
		size += domainSize[i];
	}

	size_t offset = 0;
	std::vector<std::vector<eslocal> > points(domains);
	os << "CELLS " << domains << " " << size << "\n";
	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t d = 0; d < nodes[i]->domains().size(); d++) {
			points[nodes[i]->domains()[d]].push_back(offset++);
		}
	}

	for (size_t d = 0; d < domains; d++) {
		os << points[d].size();
		for (size_t i = 0; i < points[d].size(); i++) {
			os << " " << points[d][i];
		}
		os << "\n";
	}

	os << "\n";

	os << "CELL_TYPES " << domains << "\n";
	for (size_t d = 0; d < domains; d++) {
		os << "2\n";
	}
	os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	os << "CELL_DATA " << domains << "\n";
	os << "SCALARS decomposition int 1\n";
	os << "LOOKUP_TABLE decomposition\n";
	for (size_t d = 0; d < domains; d++) {
		os << d << "\n";
	}
	os << "\n";
}

static void elements(std::ofstream &os, const espreso::Mesh &mesh)
{
	std::vector<espreso::Element*> elements(mesh.elements());
	elements.insert(elements.end(), mesh.faces().begin(), mesh.faces().end());
	elements.insert(elements.end(), mesh.edges().begin(), mesh.edges().end());

	std::sort(
			elements.begin() + mesh.elements().size(),
			elements.begin() + mesh.elements().size() + mesh.faces().size(),
			[] (const espreso::Element* e1, const espreso::Element *e2) { return e1->domains() < e2->domains(); });

	std::sort(
			elements.begin() + mesh.elements().size() + mesh.faces().size(),
			elements.end(),
			[] (const espreso::Element* e1, const espreso::Element *e2) { return e1->domains() < e2->domains(); });

	size_t nSize = 0, eSize = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		nSize += elements[i]->domains().size() * (elements[i]->nodes() + 1);
		eSize += elements[i]->domains().size();
	}

	std::vector<size_t> offset = { 0 };
	for (size_t p = 1; p < mesh.parts(); p++) {
		offset.push_back(offset[p - 1] + mesh.coordinates().localSize(p - 1));
	}

	// ELEMENTS
	os << "CELLS " << eSize << " " << nSize << "\n";
	for (size_t i = 0; i < elements.size(); i++) {
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			os << elements[i]->nodes();
			for (size_t j = 0; j < elements[i]->nodes(); j++) {
				os << " " << mesh.coordinates().localIndex(elements[i]->node(j), elements[i]->domains()[d]) + offset[elements[i]->domains()[d]];
			}
			os << "\n";
		}
	}
	os << "\n";

	// ELEMENTS TYPES
	os << "CELL_TYPES " << eSize << "\n";
	for (size_t i = 0; i < elements.size(); i++) {
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			os << elements[i]->vtkCode() << "\n";
		}
	}
	os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	os << "CELL_DATA " << eSize << "\n";
	os << "SCALARS decomposition int 1\n";
	os << "LOOKUP_TABLE decomposition\n";
	size_t part = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		if (i && elements[i]->domains() != elements[i - 1]->domains()) {
			part++;
		}
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			os << part << "\n";
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

VTK::VTK(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster): Store(mesh, path, shrinkSubdomain, shringCluster)
{
	switch (config::output::OUTPUT_FORMAT) {
	case config::output::OUTPUT_FORMATAlternatives::VTK_LEGACY_FORMAT:
		break;
	default:
		ESINFO(ALWAYS) << TextColor::YELLOW << "Warning: ESPRESO not contains a library for saving generic VTK format. VTK Legacy format is used.";
	}
	computeCenters();
}

void VTK::storeGeometry(size_t timeStep)
{
	open(_os, _path, timeStep);
	head(_os);
	coordinates(_os, _mesh.coordinates(), [&] (const Point &point, size_t part) { return shrink(point, part); });
	elements(_os, _mesh);
	_os.flush();
}

void VTK::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{
	ESINFO(GLOBAL_ERROR) << "Implement store property";
}

void VTK::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	switch (eType) {
	case ElementType::ELEMENTS:
	case ElementType::FACES:
	case ElementType::EDGES:
		ESINFO(GLOBAL_ERROR) << "Implement store values";
	case ElementType::NODES:
		coordinatesDisplacement(_os, values, dimension);
	}
	_os.flush();
}

void VTK::store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster)
{
	std::stringstream ss;
	ss << _path << config::env::MPIrank << ".vtk";

	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	head(os);
	coordinates(os, _mesh.coordinates(), [&] (const Point &point, size_t part) { return this->shrink(point, part); });
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

	VTK vtk(mesh, path, shrinkSubdomain, shrinkCluster);

	head(os);
	coordinates(os, mesh.coordinates(), [&] (const Point &point, size_t part) { return vtk.shrink(point, part); });
	elements(os, mesh);

	os.close();
}

void VTK::properties(const Mesh &mesh, const std::string &path, std::vector<Property> properties, double shrinkSubdomain, double shrinkCluster)
{
	std::stringstream ss;
	const std::vector<Element*> &elements = mesh.elements();
	const std::vector<eslocal> &partition = mesh.getPartition();

	ss << path << config::env::MPIrank << ".vtk";

	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);
	head(os);
	os << "DATASET POLYDATA\n";
	os << "POINTS " << 2*elements.size() << " float\n";
	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t e = partition[p]; e < partition[p + 1]; e++) {
			Point mid;
			for (size_t i = 0; i < elements[e]->nodes(); i++) {
				mid += mesh.coordinates().get(elements[e]->node(i), p);
			}
			mid /= elements[e]->nodes();
			os<<mid.x<<" "<<mid.y<<" "<<mid.z<<"\n";

			const std::vector<Evaluator*> &ux = elements[e]->settings(properties[0]);
			const std::vector<Evaluator*> &uy = elements[e]->settings(properties[1]);			
			
		}
	}

	os<<"LINES "<<elements.size()<<" "<<elements.size()*3<<"\n";
	for(int i=0;i<elements.size();i++){
	  os<<2<<" "<<i<<" "<<elements.size()+i<<"\n";
	}

	os<<"\nPOINT_DATA "<<elements.size()*2<<"\nVECTORS Vectors float\n";
	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t e = partition[p]; e < partition[p + 1]; e++) {
			Point mid;
			for (size_t i = 0; i < elements[e]->nodes(); i++) {
				mid += mesh.coordinates().get(elements[e]->node(i), p);
			}
			mid /= elements[e]->nodes();

			const std::vector<Evaluator*> &ux = elements[e]->settings(properties[0]);
			const std::vector<Evaluator*> &uy = elements[e]->settings(properties[1]);

			// TODO: change
			double x = 0, y = 0, z = 0;
//			double x=ux.back()->evaluate(mid.x,mid.y,mid.z)/elements[e]->nodes();
//			double y=uy.back()->evaluate(mid.x,mid.y,mid.z)/elements[e]->nodes();
			os<<x<<" "<<y<<" "<<z<<"\n";			
			
		}
	}
	//head(os);
	//coordinates(os, mesh.coordinates(),shrinkSubdomain, shrinkCluster);
	//elements(os, mesh);

	// TODO:

	os.close();
}

void VTK::fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shrinkCluster)
{
	std::vector<Element*> fixPoints;
	for (size_t p = 0; p < mesh.parts(); p++) {
		fixPoints.insert(fixPoints.end(), mesh.fixPoints(p).begin(), mesh.fixPoints(p).end());
	}

	std::sort(fixPoints.begin(), fixPoints.end());
	Esutils::removeDuplicity(fixPoints);

	std::stringstream ss;
	ss << path << config::env::MPIrank << ".vtk";
	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	VTK vtk(mesh, path, shrinkSubdomain, shrinkCluster);

	head(os);
	coordinates(os, mesh.coordinates(), fixPoints, [&] (const Point &point, size_t part) { return vtk.shrink(point, part); });
	nodes(os, fixPoints, mesh.parts());

	os.close();
}

void VTK::corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shrinkCluster)
{
	std::stringstream ss;
	ss << path << config::env::MPIrank << ".vtk";
	std::ofstream os;
	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);

	VTK vtk(mesh, path, shrinkSubdomain, shrinkCluster);

	head(os);
	coordinates(os, mesh.coordinates(), mesh.corners(), [&] (const Point &point, size_t part) { return vtk.shrink(point, part); });
	nodes(os, mesh.corners(), mesh.parts());

	os.close();
}

void VTK::gluing(const Mesh &mesh, const EqualityConstraints &constraints, const std::string &path, size_t dofs, double shrinkSubdomain, double shrinkCluster)
{

}

