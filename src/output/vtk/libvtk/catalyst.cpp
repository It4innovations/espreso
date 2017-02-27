
#include <unistd.h>

#include "../catalyst.h"

#include "../../../mesh/elements/element.h"
#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/structures/coordinates.h"

#include "../../../configuration/output.h"

#include "vtkNew.h"
#include "vtkUnstructuredGrid.h"

#include "vtkCPProcessor.h"
#include "vtkCPDataDescription.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkCPInputDataDescription.h"

using namespace espreso::store;

static void setSettings(int s, std::string name){
	vtkSmartPointer<vtkIntArray> scale =  vtkSmartPointer<vtkIntArray>::New();
	scale->SetNumberOfComponents(1);
	scale->SetName("Scale");
	scale->InsertNextValue(s);
	vtkSmartPointer<vtkStringArray> label =  vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("Label");
	label->InsertNextValue(name);
	vtkFieldData* data=vtkFieldData::New();
	data->AddArray(scale);
	data->AddArray(label);

	if (processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
		dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
		processor->CoProcess(dataDescription.GetPointer());
	}
	dataDescription->SetUserData(data);
	dataDescription->AddInput("input");
	dataDescription->SetTimeData(xx++, xx++);
	dataDescription->ForceOutputOn();
}

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: VTK(output, mesh, path), timeStep(0)
{
	processor = vtkCPProcessor::New();
	dataDescription = vtkCPDataDescription::New();

	processor->Initialize();

	vtkNew<vtkCPPythonScriptPipeline> pipeline;
	pipeline->Initialize("src/output/vtk/catalyst.py");
	processor->AddPipeline(pipeline.GetPointer());

	dataDescription->AddInput("input");
	dataDescription->ForceOutputOn();
}

Catalyst::~Catalyst()
{
	processor->Finalize();
	processor->Delete();
	dataDescription->Delete();
}

void Catalyst::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	VTK::storeValues(name, dimension, values, eType);

	dataDescription->SetTimeData(timeStep++, timeStep++);
	dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
	processor->CoProcess(dataDescription);
	sleep(_output.sleep);
}

void Catalyst::finalize()
{

}
