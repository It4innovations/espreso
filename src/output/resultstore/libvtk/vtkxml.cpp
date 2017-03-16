
#include <fstream>

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLWriter.h"

#include "../../../assembler/solution.h"
#include "../../../mesh/structures/elementtypes.h"

#include "../vtkxmlascii.h"

using namespace espreso::output;


VTKXML::VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path), _VTKGrid(NULL), _writer(NULL), _os(NULL)
{
	preprocessing();
}

VTKXML::~VTKXML()
{

}

void VTKXML::initWriter(const std::string &name, size_t points, size_t cells)
{
	_VTKGrid = vtkUnstructuredGrid::New();
	_writer->SetFileName((name + ".vtu").c_str());
	_writer->SetInputData(_VTKGrid);
}

void VTKXML::addMesh(const RegionInfo *regionInfo)
{
	// TODO: avoid copying
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetDataTypeToDouble();
	for (size_t i = 0; i < regionInfo->coordinates.size(); i += 3) {
		points->InsertNextPoint(regionInfo->coordinates[i + 0], regionInfo->coordinates[i + 1], regionInfo->coordinates[i + 2]);
	}
	_VTKGrid->SetPoints(points);

	_VTKGrid->Allocate(static_cast<vtkIdType>(regionInfo->elements.size()));

	std::vector<vtkIdType> nodes(20);
	for (size_t i = 0, p = 0; i < regionInfo->elementsTypes.size(); p = regionInfo->elementsNodes[i++]) {
		nodes.clear();
		nodes.insert(nodes.end(), regionInfo->elements.begin() + p, regionInfo->elements.begin() + regionInfo->elementsNodes[i]);
		_VTKGrid->InsertNextCell(regionInfo->elementsTypes[i], regionInfo->elementsNodes[i] - p, nodes.data());
	}
}


void VTKXML::addData(const DataArrays &data, const std::vector<Solution*> &solution)
{
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




