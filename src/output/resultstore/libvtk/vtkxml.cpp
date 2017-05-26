
#include <fstream>

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLWriter.h"

#include "../../regiondata.h"

#include "../../../assembler/solution.h"
#include "../../../mesh/structures/elementtypes.h"

#include "../vtkxmlascii.h"

using namespace espreso::output;


VTKXML::VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path), _VTKGrid(NULL), _writer(NULL), _os(NULL)
{

}

VTKXML::~VTKXML()
{

}

std::string VTKXML::initWriter(const std::string &name, size_t points, size_t cells)
{
	_VTKGrid = vtkUnstructuredGrid::New();
	_writer->SetFileName((name + ".vtu").c_str());
	_writer->SetInputData(_VTKGrid);
	return name + ".vtu";
}

void VTKXML::addMesh(const RegionData &regionData)
{
	// TODO: avoid copying
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetDataTypeToDouble();
	for (size_t i = 0; i < regionData.coordinates.size(); i += 3) {
		points->InsertNextPoint(regionData.coordinates[i + 0], regionData.coordinates[i + 1], regionData.coordinates[i + 2]);
	}
	_VTKGrid->SetPoints(points);

	_VTKGrid->Allocate(static_cast<vtkIdType>(regionData.elements.size()));

	std::vector<vtkIdType> nodes(20);
	for (size_t i = 0, p = 0; i < regionData.elementsTypes.size(); p = regionData.elementsNodes[i++]) {
		nodes.clear();
		nodes.insert(nodes.end(), regionData.elements.begin() + p, regionData.elements.begin() + regionData.elementsNodes[i]);
		_VTKGrid->InsertNextCell(regionData.elementsTypes[i], regionData.elementsNodes[i] - p, nodes.data());
	}
}


void VTKXML::addData(const RegionData &regionData)
{
	for (auto it = regionData.data.pointDataInteger.begin(); it != regionData.data.pointDataInteger.end(); ++it) {
		storePointData(it->first, it->second.first, *it->second.second);
	}
	for (auto it = regionData.data.pointDataDouble.begin(); it != regionData.data.pointDataDouble.end(); ++it) {
		storePointData(it->first, it->second.first, *it->second.second);
	}

	for (auto it = regionData.data.elementDataInteger.begin(); it != regionData.data.elementDataInteger.end(); ++it) {
		storeCellData(it->first, it->second.first, *it->second.second);
	}
	for (auto it = regionData.data.elementDataDouble.begin(); it != regionData.data.elementDataDouble.end(); ++it) {
		storeCellData(it->first, it->second.first, *it->second.second);
	}
}

void VTKXML::finalizeWriter()
{
	_writer->Write();
	_VTKGrid->Delete();
	for (size_t i = 0; i < _VTKDataArrays.size(); i++) {
		delete[] _VTKDataArrays[i];
	}
	_VTKDataArrays.clear();
}




