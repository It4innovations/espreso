
#include "../paraview.h"

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"


vtkCPProcessor* processor;
vtkSmartPointer<vtkUnstructuredGrid> VTKGrid;
vtkNew<vtkCPDataDescription> dataDescription;
double xx = 1;

using namespace espreso::output;

Paraview::Paraview(const Mesh &mesh, const std::string &path): Store(mesh, path)
{
	processor = vtkCPProcessor::New();
	processor->Initialize();

	vtkNew<vtkCPPythonScriptPipeline> pipeline;
	pipeline->Initialize("catalyst_pipeline_cube.py");
	processor->AddPipeline(pipeline.GetPointer());

	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<eslocal> &_partPtrs = _mesh.getPartition();

	VTKGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	size_t n_nodsClust = 0;
	for (size_t iEl = 0; iEl < elements.size(); iEl++) {
		n_nodsClust += elements[iEl]->nodes();
	}
	size_t cnt = 0, n_points = 0;
	for (size_t d = 0; d < _mesh.parts(); d++) {
		n_points += _mesh.coordinates().localSize(d);
	}

	//Points
	int counter = 0;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t d = 0; d < _mesh.parts(); d++) {
		for (size_t i = 0; i < _mesh.coordinates().localSize(d); i++) {
			Point xyz = _mesh.coordinates().get(i, d);
			points->InsertNextPoint(xyz.x, xyz.y, xyz.z);
			counter++;
		}
	}
	VTKGrid->SetPoints(points);

	VTKGrid->Allocate(static_cast<vtkIdType>(n_nodsClust));
	vtkIdType tmp[100]; //max number of  node

	//Cells
	size_t i = 0;
	cnt = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal ii = 0; ii < _partPtrs[part + 1] - _partPtrs[part]; ii++) {
			for (size_t j = 0; j < _mesh.elements()[i]->nodes(); j++) {
				tmp[j] = _mesh.coordinates().localIndex(elements[i]->node(j), part) + cnt;
			}
			int code = _mesh.elements()[i]->vtkCode();
			VTKGrid->InsertNextCell(code, elements[i]->nodes(), &tmp[0]);
			i++;
		}
		cnt += _mesh.coordinates().localSize(part);
	}

	float *decomposition_array = new float[_mesh.elements().size()];

	//Decomposition
	counter = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			float part_redefine = part;
			decomposition_array[counter] = part_redefine;
			counter++;
		}
	}

	vtkNew<vtkFloatArray> decomposition;
	decomposition->SetName("decomposition");
	decomposition->SetNumberOfComponents(1);
	decomposition->SetArray(decomposition_array, static_cast<vtkIdType>(elements.size()), 0);
	VTKGrid->GetCellData()->AddArray(decomposition.GetPointer());

	dataDescription->AddInput("input");
	dataDescription->SetTimeData(0, 0);
	dataDescription->ForceOutputOn();
}

void Paraview::store(std::vector<std::vector<double> > &displasment, double shrinkSubdomain, double shrinkCluster)
{
	unsigned int mycounter = 0;
	for (size_t i = 0; i < displasment.size(); i++) {
		mycounter += displasment[i].size();
	}
	double displacement_array[mycounter];
	unsigned int size = mycounter;

	mycounter = 0;
	for (size_t i = 0; i < displasment.size(); i++) {
		//note prim_solution.size = number of subdomains. all elements in prim_solution[0]
		for (size_t j = 0; j < (displasment[i].size() / 3); j++) {
			displacement_array[3 * mycounter + 0] = displasment[i][3 * j + 0];
			displacement_array[3 * mycounter + 1] = displasment[i][3 * j + 1];
			displacement_array[3 * mycounter + 2] = displasment[i][3 * j + 2];

			mycounter++;

		}
	}

	if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0) {
		// displacement array
		vtkNew<vtkDoubleArray> displacement;
		displacement->SetName("displacement");
		displacement->SetNumberOfComponents(3);
		displacement->SetNumberOfTuples(static_cast<vtkIdType>(mycounter));
		VTKGrid->GetPointData()->AddArray(displacement.GetPointer());
	}
	vtkDoubleArray* displacement = vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray("displacement"));
	double* displacementData = displacement_array;
	vtkIdType numTuples = displacement->GetNumberOfTuples();
	for (vtkIdType i = 0, counter = 0; i < numTuples; i++, counter++) {
		double values[3] = {
				displacementData[i * 3],
				displacementData[i * 3 + 1],
				xx * displacementData[i * 3 + 2]
		};
		xx += 0.1;
		displacement->SetTypedTuple(counter, values);
	}

	if (processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
		dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
		processor->CoProcess(dataDescription.GetPointer());
	}
	sleep(1);
}

void Paraview::storeGeometry(size_t timeStep = -1)
{

}

void Paraview::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{

}

void Paraview::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{

}

void Paraview::finalize()
{

}
