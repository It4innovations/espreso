
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
