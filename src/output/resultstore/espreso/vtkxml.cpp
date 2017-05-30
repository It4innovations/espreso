
#include <fstream>
#include <algorithm>

#include "../../regiondata.h"

#include "../../../configuration/environment.h"
#include "../../../assembler/solution.h"
#include "../../../mesh/structures/elementtypes.h"

#include "../vtkxml.h"

using namespace espreso::output;

VTKXML::VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path, MeshInfo::InfoMode mode)
: ResultStore(output, mesh, path, mode), _VTKGrid(NULL), _writer(NULL), _os(NULL)
{

}

VTKXML::~VTKXML()
{

}

std::string VTKXML::initWriter(const std::string &name, size_t points, size_t cells)
{
	_os = new std::ofstream();
	_os->open((name + ".vtu").c_str(), std::ios::out | std::ios::trunc);
	(*_os) << "<?xml version=\"1.0\"?>\n";
	(*_os) << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	(*_os) << "<UnstructuredGrid>\n";
	(*_os) << "<Piece NumberOfPoints=\"" << points << "\" NumberOfCells=\"" << cells << "\">\n";
	return name + ".vtu";
}

void VTKXML::addMesh(const RegionData &regionData)
{
	(*_os) << "  <Points>\n";
	storePointData("Points", 3, regionData.coordinates);
	(*_os) << "  </Points>\n";

	(*_os) << "  <Cells>\n";
	storeCellData("connectivity", 1, regionData.elements);
	storeCellData("offsets", 1, regionData.elementsNodes);
	storeCellData("types", 1, regionData.elementsTypes);
	(*_os) << "  </Cells>\n";
}

void VTKXML::addData(const RegionData &regionData)
{
	if (regionData.pointDataNames().size()) {
		(*_os) << "  <PointData Scalars=\"" << regionData.pointDataNames()[0] << "\">\n";
	} else {
		(*_os) << "  <PointData>\n";
	}
	for (auto it = regionData.data.pointDataInteger.begin(); it != regionData.data.pointDataInteger.end(); ++it) {
		storePointData(it->first, it->second.first, *it->second.second);
	}

	for (auto it = regionData.data.pointDataDouble.begin(); it != regionData.data.pointDataDouble.end(); ++it) {
		storePointData(it->first, it->second.first, *it->second.second);
	}

	(*_os) << "  </PointData>\n";

	(*_os) << "  <CellData>\n";
	for (auto it = regionData.data.elementDataInteger.begin(); it != regionData.data.elementDataInteger.end(); ++it) {
		storeCellData(it->first, it->second.first, *it->second.second);
	}

	for (auto it = regionData.data.elementDataDouble.begin(); it != regionData.data.elementDataDouble.end(); ++it) {
		storeCellData(it->first, it->second.first, *it->second.second);
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



