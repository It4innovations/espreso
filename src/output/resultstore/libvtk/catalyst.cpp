
#include <unistd.h>

#include "../catalyst.h"
#include "../../../assembler/step.h"
#include "../../../assembler/solution.h"
#include "../../../configuration/output.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

#include "vtkCPProcessor.h"
#include "vtkCPDataDescription.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkCPInputDataDescription.h"
#include "vtkFieldData.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"

using namespace espreso::output;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: VTKXMLASCII(output, mesh, path)
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
	addMesh(_coordinates, _elementsTypes, _elementsNodes, _elements);
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

void Catalyst::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{
	if (step.step == 0 && step.substep == 0 && step.iteration == 0) {
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
	}

	addData(_coordinates.size() / 3, _elementsTypes.size(), solution);
	_dataDescription->SetTimeData(step.iteration, step.iteration);
	_dataDescription->GetInputDescriptionByName("input")->SetGrid(_VTKGrid);
	_processor->CoProcess(_dataDescription);
	sleep(_configuration.sleep);
}



