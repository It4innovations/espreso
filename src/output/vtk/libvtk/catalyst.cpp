
#include "../catalyst.h"

#include "../../../mesh/elements/element.h"
#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/structures/coordinates.h"

#include "vtkNew.h"

#include "vtkDoubleArray.h"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkUnstructuredGrid.h"

#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"

using namespace espreso::store;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: ResultStore(output, mesh, path)
{
	processor = vtkCPProcessor::New();
	VTKGrid = vtkUnstructuredGrid::New();
	dataDescription = vtkCPDataDescription::New();

	processor->Initialize();

	vtkNew<vtkCPPythonScriptPipeline> pipeline;
	pipeline->Initialize("plane.py");
	processor->AddPipeline(pipeline.GetPointer());

	dataDescription->AddInput("input");
	dataDescription->SetTimeData(0, 0);
	dataDescription->ForceOutputOn();
}

Catalyst::~Catalyst()
{
	processor->Finalize();
	processor->Delete();
	VTKGrid->Delete();
	dataDescription->Delete();
	for (size_t i = 0; i < VTKDataArrays.size(); i++) {
		delete[] VTKDataArrays[i];
	}
}

void Catalyst::storeGeometry(size_t timeStep)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t d = 0; d < _mesh.parts(); d++) {
		for (size_t i = 0; i < _mesh.coordinates().localSize(d); i++) {
			const espreso::Point &p = _mesh.coordinates().get(i, d);
			points->InsertNextPoint(p.x, p.y, p.z);
		}
	}
	VTKGrid->SetPoints(points);

	size_t nSize = 0;
	for (size_t i = 0; i < _mesh.elements().size(); i++) {
		nSize += _mesh.elements()[i]->domains().size() * _mesh.elements()[i]->nodes();
	}

	std::vector<size_t> offset = { 0 };
	for (size_t p = 1; p < _mesh.parts(); p++) {
		offset.push_back(offset[p - 1] + _mesh.coordinates().localSize(p - 1));
	}

	VTKGrid->Allocate(static_cast<vtkIdType>(nSize));

	std::vector<vtkIdType> nodes(20);
	for (size_t i = 0, c = 0; i < _mesh.elements().size(); i++) {
		for (size_t d = 0; d < _mesh.elements()[i]->domains().size(); d++) {
			nodes.clear();
			for (size_t n = 0; n < _mesh.elements()[i]->nodes(); n++) {
				nodes.push_back(_mesh.coordinates().localIndex(_mesh.elements()[i]->node(n), _mesh.elements()[i]->domains()[d]) + offset[_mesh.elements()[i]->domains()[d]]);
			}
			VTKGrid->InsertNextCell(_mesh.elements()[i]->vtkCode(), _mesh.elements()[i]->nodes(), nodes.data());
		}
	}
}

void Catalyst::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{
	ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Catalyst skips storing properties.";
}

void Catalyst::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	size_t size = 0;
	for (size_t i = 0; i < values.size(); i++) {
		size += values[i].size();
	}

	double *data = new double[size];
	for (size_t i = 0, offset = 0; i < values.size(); offset += values[i++].size()) {
		memcpy(data + offset, values[i].data(), values[i].size() * sizeof(double));
	}

	vtkNew<vtkDoubleArray> vtkArray;
	vtkArray->SetName(name.c_str());
	vtkArray->SetNumberOfComponents(dimension);
	vtkArray->SetArray(data, static_cast<vtkIdType>(size), 1);

	switch (eType) {
	case espreso::store::ElementType::NODES:
		VTKGrid->GetPointData()->AddArray(vtkArray.GetPointer());
		if (VTKGrid->GetPointData()->GetNumberOfArrays() == 1) {
			VTKGrid->GetPointData()->SetActiveScalars(name.c_str());
		}
		break;
	case espreso::store::ElementType::ELEMENTS:
		VTKGrid->GetCellData()->AddArray(vtkArray.GetPointer());
		if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0) {
			VTKGrid->GetCellData()->SetActiveScalars(name.c_str());
		}
		break;
	}
	VTKDataArrays.push_back(data);

	dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
	processor->CoProcess(dataDescription);
}

void Catalyst::finalize()
{

}
