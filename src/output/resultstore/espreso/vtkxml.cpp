
#include <fstream>
#include <algorithm>

#include "../../../configuration/environment.h"
#include "../../../assembler/solution.h"
#include "../../../mesh/structures/elementtypes.h"

#include "../vtkxml.h"

using namespace espreso::output;

VTKXML::VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path), _VTKGrid(NULL), _os(NULL)
{
	preprocessing();
	size_t offset = 0;
	std::for_each(_elementsNodes.begin(), _elementsNodes.end(), [&] (eslocal &nodes) { nodes = offset += nodes; });
}

void VTKXML::regionPreprocessing(const espreso::Region &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	ResultStore::regionPreprocessing(region, coordinates, elementsTypes, elementsNodes, elements);
	size_t offset = 0;
	std::for_each(elementsNodes.begin(), elementsNodes.end(), [&] (eslocal &nodes) { nodes = offset += nodes; });
}


VTKXML::~VTKXML()
{

}

void VTKXML::finalize()
{

}

void VTKXML::initWriter(const std::string &name, size_t points, size_t cells)
{
	_os = new std::ofstream();
	_os->open((name + ".vtu").c_str(), std::ios::out | std::ios::trunc);
	(*_os) << "<?xml version=\"1.0\"?>\n";
	(*_os) << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	(*_os) << "<UnstructuredGrid>\n";
	(*_os) << "<Piece NumberOfPoints=\"" << points << "\" NumberOfCells=\"" << cells << "\">\n";
}

void VTKXML::addMesh(std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	(*_os) << "  <Points>\n";
	(*_os) << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << format() << "\">\n";
	(*_os) << "      "; store(*_os, coordinates); *_os << "\n";
	(*_os) << "    </DataArray>\n";
	(*_os) << "  </Points>\n";

	(*_os) << "  <Cells>\n";
	(*_os) << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"" << format() << "\">\n";
	(*_os) << "      "; store(*_os, elements); *_os << "\n";
	(*_os) << "    </DataArray>\n";
	(*_os) << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"" << format() << "\">\n";
	(*_os) << "      "; store(*_os, elementsNodes); *_os << "\n";
	(*_os) << "    </DataArray>\n";
	(*_os) << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"" << format() << "\">\n";
	(*_os) << "      "; store(*_os, elementsTypes); *_os << "\n";
	(*_os) << "    </DataArray>\n";
	(*_os) << "  </Cells>\n";
}

void VTKXML::addData(size_t points, size_t cells, const DataArrays &data)
{
	(*_os) << "  <PointData>\n";
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		(*_os) << "    <DataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() / points << "\">\n";
		(*_os) << "      "; store((*_os), *it->second); (*_os) << "\n";
		(*_os) << "    </DataArray>\n";
	}

	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		(*_os) << "    <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() / points << "\">\n";
		(*_os) << "      "; store((*_os), *it->second); (*_os) << "\n";
		(*_os) << "    </DataArray>\n";
	}
	(*_os) << "  </PointData>\n";

	(*_os) << "  <CellData>\n";
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		(*_os) << "    <DataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() / cells << "\">\n";
		(*_os) << "      "; store((*_os), *it->second); (*_os) << "\n";
		(*_os) << "    </DataArray>\n";
	}

	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		(*_os) << "    <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() /cells << "\">\n";
		(*_os) << "      "; store((*_os), *it->second); (*_os) << "\n";
		(*_os) << "    </DataArray>\n";
	}
	(*_os) << "  </CellData>\n";
}

void VTKXML::addData(size_t points, size_t cells, const std::vector<Solution*> &solution)
{
	(*_os) << "  <PointData>\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::NODES) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= points;
			(*_os) << "    <DataArray type=\"Float64\" Name=\"" << solution[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << components  << "\">\n";
			(*_os) << "      ";
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { store((*_os), part); });;
			(*_os) << "\n";
			(*_os) << "    </DataArray>\n";
		}
	}
	(*_os) << "  </PointData>\n";

	(*_os) << "  <CellData>\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::ELEMENTS) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= cells;
			(*_os) << "    <DataArray type=\"Float64\" Name=\"" << solution[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << components << "\">\n";
			(*_os) << "      ";
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { store((*_os), part); });;
			(*_os) << "\n";
			(*_os) << "    </DataArray>\n";
		}
	}
	(*_os) << "  </CellData>\n";
}

void VTKXML::finalizeWriter()
{
	(*_os) << "</Piece>\n";
	(*_os) << "</UnstructuredGrid>\n";
	(*_os) << "</VTKFile>\n";
	_os->close();
	delete _os;
}



