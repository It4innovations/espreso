
#include <unistd.h>

#include "../catalyst.h"
#include "../../../assembler/step.h"
#include "../../../assembler/solution.h"
#include "../../../configuration/output.h"

#include "../../regiondata.h"
#include "../../meshinfo.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

#include "vtkCPProcessor.h"
#include "vtkCPDataDescription.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkCPInputDataDescription.h"
#include "vtkFieldData.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"

using namespace espreso;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode)
: VTKXMLASCII(output, mesh, mode)
{
	_processor = vtkCPProcessor::New();
	_dataDescription = vtkCPDataDescription::New();
	_fieldData = vtkFieldData::New();

	_processor->Initialize();

	vtkSmartPointer<vtkCPPythonScriptPipeline> pipeline = vtkSmartPointer<vtkCPPythonScriptPipeline>::New();
	pipeline->Initialize("src/output/resultstore/catalyst.py");
	_processor->AddPipeline(pipeline.GetPointer());

	_dataDescription->AddInput("input");
	_dataDescription->ForceOutputOn();

	_VTKGrid = vtkUnstructuredGrid::New();
	addMesh(_meshInfo->region(0)); // TODO: show more bodies via Catalyst
}

Catalyst::~Catalyst()
{
	_VTKGrid->Delete();
	_processor->Finalize();
	_processor->Delete();
	_dataDescription->Delete();
	_fieldData->Delete();
}

void Catalyst::storeSettings(const Step &step)
{
	// only solution can be stored by catalyst
}

void Catalyst::storeSettings(size_t steps)
{
	// only solution can be stored by catalyst
}

void Catalyst::storeSettings(const std::vector<size_t> &steps)
{
	// only solution can be stored by catalyst
}

void Catalyst::storeFETIData(const Step &step, const Instance &instance)
{
	// only solution can be stored by catalyst
}

void Catalyst::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{
	vtkSmartPointer<vtkIntArray> scale = vtkSmartPointer<vtkIntArray>::New();
	scale->SetNumberOfComponents(1);
	scale->SetName("scale");
	scale->InsertNextValue(2);
	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");
	label->InsertNextValue(solution.front()->name);

	_fieldData->AddArray(scale);
	_fieldData->AddArray(label);
	_dataDescription->SetUserData(_fieldData);

	_meshInfo->addSolution(solution);
	addData(_meshInfo->region(0));

	_dataDescription->SetTimeData(step.substep, step.substep);
	_dataDescription->GetInputDescriptionByName("input")->SetGrid(_VTKGrid);
	_processor->CoProcess(_dataDescription);
	sleep(_configuration.sleep);
}



