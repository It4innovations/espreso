
#include "../../../old_output/vtk/catalyst.h"

#include <unistd.h>

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
#include "vtkFieldData.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"

using namespace espreso::store;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: VTK(output, mesh, path), timeStep(0)
{
	processor = vtkCPProcessor::New();
	dataDescription = vtkCPDataDescription::New();
	fieldData = vtkFieldData::New();

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
	fieldData->Delete();
}

void Catalyst::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	VTK::storeValues("results", dimension, values, eType);

	if (timeStep == 0) {
		vtkSmartPointer<vtkIntArray> scale = vtkSmartPointer<vtkIntArray>::New();
		scale->SetNumberOfComponents(1);
		scale->SetName("scale");
		scale->InsertNextValue(2);
		vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
		label->SetNumberOfComponents(1);
		label->SetName("label");
		label->InsertNextValue(name);

		fieldData->AddArray(scale);
		fieldData->AddArray(label);
		dataDescription->SetUserData(fieldData);
	}

	dataDescription->SetTimeData(timeStep++, timeStep++);
	dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
	processor->CoProcess(dataDescription);
	sleep(_output.sleep);
}

void Catalyst::finalize()
{

}
