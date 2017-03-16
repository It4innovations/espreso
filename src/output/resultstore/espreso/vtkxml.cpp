
#include <fstream>
#include <algorithm>

#include "../../../configuration/environment.h"
#include "../../../assembler/solution.h"
#include "../../../mesh/structures/elementtypes.h"

#include "../vtkxml.h"

using namespace espreso::output;

VTKXML::VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path), _VTKGrid(NULL), _writer(NULL), _os(NULL)
{

}

VTKXML::~VTKXML()
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

void VTKXML::addMesh(const RegionInfo *regionInfo)
{
	(*_os) << "  <Points>\n";
	storePointData("Points", 3, regionInfo->coordinates);
	(*_os) << "  </Points>\n";

	(*_os) << "  <Cells>\n";
	storeCellData("connectivity", 1, regionInfo->elements);
	storeCellData("offsets", 1, regionInfo->elementsNodes);
	storeCellData("types", 1, regionInfo->elementsTypes);
	(*_os) << "  </Cells>\n";
}

void VTKXML::addData(const DataArrays &data, const std::vector<Solution*> &solution)
{
	(*_os) << "  <PointData>\n";
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		storePointData(it->first, it->second.first, *it->second.second);
	}

	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		storePointData(it->first, it->second.first, *it->second.second);
	}

	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::NODES) {
			storePointData(solution[i]->name, solution[i]->properties, solution[i]->data);
		}
	}
	(*_os) << "  </PointData>\n";

	(*_os) << "  <CellData>\n";
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		storeCellData(it->first, it->second.first, *it->second.second);
	}

	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		storeCellData(it->first, it->second.first, *it->second.second);
	}
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::ELEMENTS) {
			storeCellData(solution[i]->name, solution[i]->properties, solution[i]->data);
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



