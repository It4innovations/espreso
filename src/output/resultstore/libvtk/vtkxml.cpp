
#include <fstream>
#include <algorithm>

#include "vtkNew.h"
#include "vtkSmartPointer.h"

#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkDoubleArray.h"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkUnstructuredGrid.h"

#include "vtkGeometryFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkTriangleFilter.h"
#include "vtkDecimatePro.h"
#include "vtkAppendFilter.h"

#include "vtkXMLUnstructuredGridWriter.h"

#include "../../../configuration/environment.h"
#include "../../../assembler/solution.h"
#include "../../../mesh/structures/elementtypes.h"

#include "../vtkxmlascii.h"

using namespace espreso::output;

VTKXML::VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path)
{
	preprocessing();
	_VTKGrid = vtkUnstructuredGrid::New();
}

VTKXML::~VTKXML()
{
	_VTKGrid->Delete();
}

void VTKXML::finalize()
{

}

void VTKXML::regionPreprocessing(const espreso::Region &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	ResultStore::regionPreprocessing(region, coordinates, elementsTypes, elementsNodes, elements);
}

void VTKXML::storeMesh(std::ofstream &os, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	// TODO: avoid copying
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t i = 0; i < coordinates.size(); i += 3) {
		points->InsertNextPoint(coordinates[i + 0], coordinates[i + 1], coordinates[i + 2]);
	}
	_VTKGrid->SetPoints(points);

	_VTKGrid->Allocate(static_cast<vtkIdType>(elements.size()));

	std::vector<vtkIdType> nodes(20);
	for (size_t i = 0, p = 0; i < elementsTypes.size(); p += elementsNodes[i++]) {
		nodes.clear();
		nodes.insert(nodes.end(), elements.begin() + p, elements.begin() + p + elementsNodes[i]);
		_VTKGrid->InsertNextCell(elementsTypes[i], elementsNodes[i], nodes.data());
	}
}

template <typename TVTKType, typename TType>
static void addPointArray(vtkUnstructuredGrid *VTKGrid, const std::string &name, std::vector<TType> &data)
{
	vtkNew<TVTKType> vtkArray;
	vtkArray->SetName(name.c_str());
	vtkArray->SetNumberOfComponents(1);
	vtkArray->SetArray(data.data(), static_cast<vtkIdType>(data.size()), 1);

	VTKGrid->GetPointData()->AddArray(vtkArray.GetPointer());
	if (VTKGrid->GetPointData()->GetNumberOfArrays() == 1) {
		VTKGrid->GetPointData()->SetActiveScalars(name.c_str());
	}
}

template <typename TVTKType, typename TType>
static void addCellArray(vtkUnstructuredGrid *VTKGrid, const std::string &name, std::vector<TType> &data)
{
	vtkNew<TVTKType> vtkArray;
	vtkArray->SetName(name.c_str());
	vtkArray->SetNumberOfComponents(1);
	vtkArray->SetArray(data.data(), static_cast<vtkIdType>(data.size()), 1);

	VTKGrid->GetCellData()->AddArray(vtkArray.GetPointer());
	if (VTKGrid->GetCellData()->GetNumberOfArrays() == 1) {
		VTKGrid->GetCellData()->SetActiveScalars(name.c_str());
	}
}

static void addPointArray(vtkUnstructuredGrid *VTKGrid, const std::string &name, std::vector<int> &data)
{
	addPointArray<vtkIntArray, int>(VTKGrid, name, data);
}

static void addPointArray(vtkUnstructuredGrid *VTKGrid, const std::string &name, std::vector<long> &data)
{
	addPointArray<vtkLongArray, long>(VTKGrid, name, data);
}

static void addCellArray(vtkUnstructuredGrid *VTKGrid, const std::string &name, std::vector<int> &data)
{
	addCellArray<vtkIntArray, int>(VTKGrid, name, data);
}

static void addCellArray(vtkUnstructuredGrid *VTKGrid, const std::string &name, std::vector<long> &data)
{
	addCellArray<vtkLongArray, long>(VTKGrid, name, data);
}

void VTKXML::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data)
{
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		addPointArray(_VTKGrid, it->first, *it->second);
	}
	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		addPointArray<vtkDoubleArray>(_VTKGrid, it->first, *it->second);
	}

	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		addCellArray(_VTKGrid, it->first, *it->second);
	}
	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		addCellArray<vtkDoubleArray>(_VTKGrid, it->first, *it->second);
	}

	writer->Write();
}

void VTKXML::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution)
{

}

void VTKXML::savePVTU(const std::string &root, const std::string &name, const std::vector<std::pair<std::string, size_t> > &pointDataInt, const std::vector<std::pair<std::string, size_t>> &pointDataDouble, const std::vector<std::pair<std::string, size_t>> &cellDataInt, const std::vector<std::pair<std::string, size_t>> &cellDataDouble)
{
	std::ofstream os;

	os.open(root + name + ".pvtu", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
	os << "\n";
	os << "    <PPoints>\n";
	os << "      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << format() << "\"/>\n";
	os << "    </PPoints>\n";
	os << "\n";
	os << "    <PPointData>\n";
	for (size_t i = 0; i < pointDataInt.size(); i++) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << pointDataInt[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << pointDataInt[i].second << "\"/>\n";
	}
	for (size_t i = 0; i < pointDataDouble.size(); i++) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << pointDataDouble[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << pointDataDouble[i].second << "\"/>\n";
	}
	os << "    </PPointData>\n";
	os << "\n";
	os << "    <PCellData>\n";
	for (size_t i = 0; i < cellDataInt.size(); i++) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << cellDataInt[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << cellDataInt[i].second << "\"/>\n";
	}
	for (size_t i = 0; i < cellDataDouble.size(); i++) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << cellDataDouble[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << cellDataDouble[i].second << "\"/>\n";
	}
	os << "    </PCellData>\n";
	os << "\n";
	for (int r = 0; r < environment->MPIsize; r++) {
		os << "    <Piece Source=\"" << r << "/" << name << ".vtu\"/>\n";
	}
	os << "  </PUnstructuredGrid>\n";
	os << "</VTKFile>\n";
}

void VTKXML::composeClusters(const std::string &root, const std::string &name, const DataArrays &data)
{
	std::vector<std::pair<std::string, size_t> > pointDataInt, pointDataDouble, cellDataInt, cellDataDouble;

	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		cellDataInt.push_back(std::make_pair(it->first, it->second->size() / _elementsTypes.size()));
	}
	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		cellDataDouble.push_back(std::make_pair(it->first, it->second->size() / _elementsTypes.size()));
	}
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		pointDataInt.push_back(std::make_pair(it->first, it->second->size() / (_coordinates.size() / 3)));
	}
	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		pointDataDouble.push_back(std::make_pair(it->first, it->second->size() / (_coordinates.size() / 3)));
	}
	savePVTU(root, name, pointDataInt, pointDataDouble, cellDataInt, cellDataDouble);
}

void VTKXML::composeClusters(const std::string &root, const std::string &name, const std::vector<Solution*> &solution)
{
	std::vector<std::pair<std::string, size_t> > pointData, cellData;

	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::ELEMENTS) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= _elementsTypes.size();
			cellData.push_back(std::make_pair(solution[i]->name, components));
		}
		if (solution[i]->eType == ElementType::NODES) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= (_coordinates.size() / 3);
			pointData.push_back(std::make_pair(solution[i]->name, components));
		}
	}
	savePVTU(root, name, {}, pointData, {}, cellData);
}

void VTKXML::composeRegions(const std::string &name, const std::vector<std::string> &names)
{
	std::ofstream os;

	os.open(name + ".vtm", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "<vtkMultiBlockDataSet>\n";
	for (size_t i = 0; i < names.size(); i++) {
		os << "  <DataSet index=\"" << i << "\" name=\"" << names[i] << "\" file=\"" << names[i] << ".pvtu\"> </DataSet>\n";
	}
	os << "</vtkMultiBlockDataSet>\n";
	os << "</VTKFile>\n";

}



