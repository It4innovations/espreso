
#include "../vtk.h"

#include "../../../basis/logging/logging.h"

#include "../../../config/output.h"
#include "../../../config/environment.h"

#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/elements/element.h"

#include "../../../mesh/elements/line/line2.h"

using namespace espreso::store;

static void create(std::ofstream &os, const std::string &path)
{
	std::stringstream ss;
	ss << path << espreso::environment->MPIrank;
	ss << ".vtk";

	os.open(ss.str().c_str(), std::ios::out | std::ios::trunc);
	os << "# vtk DataFile Version 4.0\n";
	os << "ESPRESO output\n";
	os << "ASCII\n";
	os << "\n";
}

VTK::VTK(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: Store(output, mesh, path), _lastData(ElementType::ELEMENTS)
{
	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:
		break;
	default:
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Warning: ESPRESO not contains a library for saving generic VTK format. VTK Legacy format is used.";
	}

	computeCenters();
}

VTK::~VTK()
{
}

static size_t coordinateSize(const espreso::Coordinates &coordinates)
{
	size_t cSize = 0;
	for (size_t p = 0; p < coordinates.parts(); p++) {
		cSize += coordinates.localToCluster(p).size();
	}
	return cSize;
}

void VTK::coordinates()
{
	create(_os, _path);
	size_t parts = _mesh.parts();
	size_t cSize = coordinateSize(_mesh.coordinates());

	_os << "DATASET UNSTRUCTURED_GRID\n";
	_os << "POINTS " << cSize << " float\n";

	for (size_t p = 0; p < parts; p++) {
		for (size_t i = 0; i < _mesh.coordinates().localToCluster(p).size(); i++) {
			_os << shrink(_mesh.coordinates().get(i, p), p) << "\n";
		}
	}
	_os << "\n";
	_os.flush();
}

void VTK::nodes(const std::vector<Element*> &nodes)
{
	create(_os, _path);
	std::vector<size_t> domainSize(_mesh.parts(), 0);
	size_t size = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t d = 0; d < nodes[i]->domains().size(); d++) {
			domainSize[nodes[i]->domains()[d]]++;
		}
	}

	for (size_t i = 0; i < _mesh.parts(); i++) {
		size += domainSize[i];
	}

	size_t offset = 0;
	std::vector<std::vector<eslocal> > points(_mesh.parts());
	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t d = 0; d < nodes[i]->domains().size(); d++) {
			points[nodes[i]->domains()[d]].push_back(offset++);
		}
	}

	_os << "DATASET UNSTRUCTURED_GRID\n";
	_os << "POINTS " << size << " float\n";
	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t d = 0; d < nodes[i]->domains().size(); d++) {
			_os << shrink(_mesh.coordinates()[nodes[i]->node(0)], nodes[i]->domains()[d]) << "\n";
		}
	}
	_os << "\n";

	_os << "CELLS " << _mesh.parts() << " " << size + _mesh.parts() << "\n";
	for (size_t d = 0; d < _mesh.parts(); d++) {
		_os << points[d].size();
		for (size_t i = 0; i < points[d].size(); i++) {
			_os << " " << points[d][i];
		}
		_os << "\n";
	}

	_os << "\n";

	_os << "CELL_TYPES " << _mesh.parts() << "\n";
	for (size_t d = 0; d < _mesh.parts(); d++) {
		_os << "2\n";
	}
	_os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	_os << "CELL_DATA " << _mesh.parts() << "\n";
	_os << "SCALARS decomposition int 1\n";
	_os << "LOOKUP_TABLE decomposition\n";
	for (size_t d = 0; d < _mesh.parts(); d++) {
		_os << d << "\n";
	}
	_os << "\n";
}

void VTK::nodes(const std::vector<std::vector<eslocal> > &nodes)
{
	create(_os, _path);
	std::vector<size_t> domainSize(_mesh.parts(), 0);
	size_t size = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		domainSize[i] = nodes[i].size();
	}

	for (size_t i = 0; i < _mesh.parts(); i++) {
		size += domainSize[i];
	}

	_os << "DATASET UNSTRUCTURED_GRID\n";
	_os << "POINTS " << size << " float\n";
	for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t i = 0; i < nodes[p].size(); i++) {
			_os << shrink(_mesh.coordinates()[nodes[p][i]], p) << "\n";
		}
	}
	_os << "\n";

	_os << "CELLS " << _mesh.parts() << " " << size + _mesh.parts() << "\n";

	for (size_t d = 0, index = 0; d < _mesh.parts(); d++) {
		_os << nodes[d].size();
		for (size_t i = 0; i < nodes[d].size(); i++) {
			_os << " " << index++;
		}
		_os << "\n";
	}

	_os << "\n";

	_os << "CELL_TYPES " << _mesh.parts() << "\n";
	for (size_t d = 0; d < _mesh.parts(); d++) {
		_os << "2\n";
	}
	_os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	_os << "CELL_DATA " << _mesh.parts() << "\n";
	_os << "SCALARS decomposition int 1\n";
	_os << "LOOKUP_TABLE decomposition\n";
	for (size_t d = 0; d < _mesh.parts(); d++) {
		_os << d << "\n";
	}
	_os << "\n";
}

void VTK::cells(ElementType eType)
{
	std::vector<espreso::Element*> elements;

	switch (eType) {
	case espreso::store::Store::ElementType::ELEMENTS:
		elements.insert(elements.end(), _mesh.elements().begin(), _mesh.elements().end());
		break;
	case espreso::store::Store::ElementType::FACES:
		elements.insert(elements.end(), _mesh.faces().begin(), _mesh.faces().end());
		std::sort(elements.begin(), elements.end(), [] (const espreso::Element* e1, const espreso::Element *e2) { return e1->domains() < e2->domains(); });
		break;
	case espreso::store::Store::ElementType::EDGES:
		elements.insert(elements.end(), _mesh.edges().begin(), _mesh.edges().end());
		std::sort(elements.begin(), elements.end(), [] (const espreso::Element* e1, const espreso::Element *e2) { return e1->domains() < e2->domains(); });
		break;
	default:
		break;
	}

	size_t nSize = 0, eSize = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		nSize += elements[i]->domains().size() * (elements[i]->nodes() + 1);
		eSize += elements[i]->domains().size();
	}

	std::vector<size_t> offset = { 0 };
	for (size_t p = 1; p < _mesh.parts(); p++) {
		offset.push_back(offset[p - 1] + _mesh.coordinates().localSize(p - 1));
	}

	// ELEMENTS
	_os << "CELLS " << eSize << " " << nSize << "\n";
	for (size_t i = 0; i < elements.size(); i++) {
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			_os << elements[i]->nodes();
			for (size_t j = 0; j < elements[i]->nodes(); j++) {
				_os << " " << _mesh.coordinates().localIndex(elements[i]->node(j), elements[i]->domains()[d]) + offset[elements[i]->domains()[d]];
			}
			_os << "\n";
		}
	}
	_os << "\n";

	// ELEMENTS TYPES
	_os << "CELL_TYPES " << eSize << "\n";
	for (size_t i = 0; i < elements.size(); i++) {
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			_os << elements[i]->vtkCode() << "\n";
		}
	}
	_os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	_os << "CELL_DATA " << eSize << "\n";
	_os << "SCALARS decomposition int 1\n";
	_os << "LOOKUP_TABLE decomposition\n";
	size_t part = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		if (i && elements[i]->domains() != elements[i - 1]->domains()) {
			part++;
		}
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			_os << part << "\n";
		}
	}
	_os << "\n";
	_os.flush();
}

void VTK::lambdas(const std::vector<std::vector<eslocal> > &nodes, std::function<Point(const Point&, size_t, size_t, bool)> shrink)
{
	create(_os, _path);

	size_t cSize = 0;
	for (size_t p = 0; p < nodes.size(); p++) {
		cSize += nodes[p].size();
	}

	_os << "DATASET UNSTRUCTURED_GRID\n";
	_os << "POINTS " << 2 * cSize << " float\n";

	for (size_t p = 0; p < nodes.size(); p++) {
		for (size_t n = 0; n < nodes[p].size(); n++) {
			Point p1 = shrink(_mesh.coordinates()[nodes[p][n]], p, n, true);
			Point p2 = shrink(_mesh.coordinates()[nodes[p][n]], p, n, false);
			_os << p1 << " " << p2 << "\n";
		}
	}

	// ELEMENTS
	_os << "\nCELLS " << cSize << " " << 3 * cSize << "\n";
	for (size_t p = 0, index = 0; p < nodes.size(); p++) {
		for (size_t n = 0; n < nodes[p].size(); n++, index++) {
			_os << "2 " << 2 * index << " " << 2 * index + 1 << "\n";
		}
	}
	_os << "\n";

	// ELEMENTS TYPES
	_os << "CELL_TYPES " << cSize << "\n";
	for (size_t p = 0; p < nodes.size(); p++) {
		for (size_t n = 0; n < nodes[p].size(); n++) {
			_os << Line2VTKCode << "\n";
		}
	}
	_os << "\n";

	// DECOMPOSITION TO SUBDOMAINS
	_os << "CELL_DATA " << cSize << "\n";
	_os << "SCALARS decomposition int 1\n";
	_os << "LOOKUP_TABLE decomposition\n";
	for (size_t i = 0; i < nodes.size(); i++) {
		for (size_t j = 0; j < nodes[i].size(); j++) {
			_os << i << "\n";
		}
	}
	_os << "\n";
	_os.flush();
}

template <typename Ttype>
static void setValueHeader(std::ofstream &os, const std::vector<std::vector<Ttype> > &values, size_t dimension, espreso::store::Store::ElementType &last, espreso::store::Store::ElementType current)
{
	if (last != current) {
		size_t size = 0;
		for (size_t p = 0; p < values.size(); p++) {
			size += values[p].size() / dimension;
		}
		switch (current) {
		case espreso::store::Store::ElementType::NODES:
			os << "POINT_DATA " << size << "\n";
			last = espreso::store::Store::ElementType::NODES;
			break;
		default:
			os << "CELL_DATA " << size << "\n";
			last = espreso::store::Store::ElementType::ELEMENTS;
		}
	}
}

template <typename Ttype>
static void storeData(std::ofstream &os, const std::vector<std::vector<Ttype> > &values, size_t DOFs)
{
	for (size_t p = 0; p < values.size(); p++) {
		for (size_t i = 0; i < values[p].size() / DOFs; i++) {
			for (size_t d = 0; d < DOFs; d++) {
				os << values[p][DOFs * i + d] << " ";
			}
			os << "\n";
		}
	}
	os.flush();
}

void VTK::data(const std::string &name, size_t dimension, const std::vector<std::vector<int> > &values, espreso::store::Store::ElementType eType)
{
	setValueHeader(_os, values, dimension, _lastData, eType);
	_os << "SCALARS " << name << " int " << dimension << "\n";
	_os << "LOOKUP_TABLE default\n";
	storeData<int>(_os, values, dimension);
}

void VTK::data(const std::string &name, size_t dimension, const std::vector<std::vector<long> > &values, espreso::store::Store::ElementType eType)
{
	setValueHeader(_os, values, dimension, _lastData, eType);
	_os << "SCALARS " << name << " long " << dimension << "\n";
	_os << "LOOKUP_TABLE default\n";
	storeData<long>(_os, values, dimension);
}

void VTK::data(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, espreso::store::Store::ElementType eType)
{
	setValueHeader(_os, values, dimension, _lastData, eType);
	_os << "SCALARS " << name << " double " << dimension << "\n";
	_os << "LOOKUP_TABLE default\n";
	storeData<double>(_os, values, dimension);
}

void VTK::finalize()
{
	_os.close();
}

