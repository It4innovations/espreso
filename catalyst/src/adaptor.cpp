#include "adaptor.h"

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

namespace {
vtkCPProcessor* Processor = NULL;
vtkUnstructuredGrid* VTKGrid;

void BuildVTKGrid(mesh::Mesh *mesh, std::vector<std::vector<eslocal> > &l2g_vec)
{
	//UNDER CONSTRUCTION
	int MPIsize = 1;
	int MPIrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << "##############################################here\n";
	std::cout << "MPIsize = " << MPIsize << "\n";
	std::cout << "MPIrank = " << MPIrank << "\n";
	std::cout << "##############################################here\n";
	MPI_Barrier(MPI_COMM_WORLD);

	const std::vector<mesh::Element*> &elements = mesh->getElements();
	const std::vector<eslocal> &_partPtrs = mesh->getPartition();

	size_t n_nodsClust = 0;
	for (size_t iEl = 0; iEl < elements.size(); iEl++) {
		n_nodsClust += elements[iEl]->size();
	}
	size_t nSubClst = l2g_vec.size();
	size_t cnt = 0, n_points = 0;
	for (size_t d = 0; d < nSubClst; d++) {
		n_points += l2g_vec[d].size();
	}
	double shrinking = 0.90;
	mesh::Coordinates &_coordinates = mesh->coordinates();
	double *coord_xyz = new double[n_points * 3];   // TODO not deleted

	int counter = 0;
	for (size_t d = 0; d < nSubClst; d++) {
		mesh::Point center;
		for (size_t c = 0; c < l2g_vec[d].size(); c++) {
			center += _coordinates[l2g_vec[d][c]];
		}
		center /= l2g_vec[d].size();

		for (size_t i = 0; i < l2g_vec[d].size(); i++) {
			mesh::Point xyz = _coordinates[l2g_vec[d][i]];
			xyz = center + (xyz - center) * shrinking;
			coord_xyz[3 * counter + 0] = xyz.x;
			coord_xyz[3 * counter + 1] = xyz.y;
			coord_xyz[3 * counter + 2] = xyz.z;
			counter++;
		}
	}

	vtkNew<vtkDoubleArray> pointArray;
	pointArray->SetNumberOfComponents(3);
	//
	size_t numpoints = n_points * 3;
	pointArray->SetArray(coord_xyz, numpoints, 1);
	vtkNew<vtkPoints> points;
	points->SetData(pointArray.GetPointer());
	VTKGrid->SetPoints(points.GetPointer());

	VTKGrid->Allocate(static_cast<vtkIdType>(n_nodsClust));
	vtkIdType tmp[MAX_NUMBER_OF_NODES];/* max number of nodes predicted to be '20' */
	//
	size_t i = 0;
	cnt = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal ii = 0; ii < _partPtrs[part + 1] - _partPtrs[part]; ii++) {
			for (size_t j = 0; j < elements[i]->size(); j++) {
				tmp[j] = elements[i]->node(j) + cnt;
			}
			VTKGrid->InsertNextCell(elements[i]->vtkCode(), elements[i]->size(), &tmp[0]);
			i++;
		}
		cnt += l2g_vec[part].size();
	}
	//
	float *decomposition_array = new float[elements.size()];

	counter = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			float part_redefine = (float) part;
			decomposition_array[counter] = part_redefine;
			counter++;
		}
	}
	//
	if (VTKGrid->GetCellData()->GetNumberOfArrays() == 0) {
		vtkNew<vtkFloatArray> decomposition;
		decomposition->SetName("decomposition");
		decomposition->SetNumberOfComponents(1);
		VTKGrid->GetCellData()->AddArray(decomposition.GetPointer());
	}

	vtkFloatArray* decomposition = vtkFloatArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("decomposition"));
	// The decomposition array is a scalar array so we can reuse
	// memory as long as we ordered the points properly.
	//    float* decompositionData = decomposition_array;
	decomposition->SetArray(decomposition_array, static_cast<vtkIdType>(elements.size()), 1);

}

void UpdateVTKAttributes(std::vector<std::vector<double> >& prim_solution)
{
	int MPIsize;
	int MPIrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	unsigned int mycounter = 0;
	//	if (MPIrank == 0 )
	//	{
	for (size_t i = 0; i < prim_solution.size(); i++) {
		mycounter += prim_solution[i].size();
	}
	//		MPI_Bcast(&mycounter, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	//	}
	//	if (MPIrank != 0)
	//	{
	//		MPI_Bcast(&mycounter, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	//	}
	//	MPI_Barrier(MPI_COMM_WORLD);
	double displacement_array[mycounter];
	unsigned int size = mycounter;

	mycounter = 0;
	for (size_t i = 0; i < prim_solution.size(); i++) {
		//note prim_solution.size = number of subdomains. all elements in prim_solution[0]
		for (size_t j = 0; j < (prim_solution[i].size() / 3); j++) {
			displacement_array[3 * mycounter + 0] = prim_solution[i][3 * j + 0];
			displacement_array[3 * mycounter + 1] = prim_solution[i][3 * j + 1];
			displacement_array[3 * mycounter + 2] = prim_solution[i][3 * j + 2];

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
				displacementData[i * 3 + 2]
		};
		displacement->SetTupleValue(counter, values);
	}
}

void BuildVTKDataStructures(
		mesh::Mesh *mesh,
		std::vector<std::vector<eslocal> > &l2g_vec,
		std::vector<std::vector<double> >& prim_solution)
{
	if (VTKGrid == NULL) {
		// The grid structure isn't changing so we only build it
		// the first time it's needed. If we needed the memory
		// we could delete it and rebuild as necessary.
		VTKGrid = vtkUnstructuredGrid::New();
		BuildVTKGrid(mesh, l2g_vec);
	}
	UpdateVTKAttributes(prim_solution);	//, decomposition_values);
}

}	// end of namespace

namespace Adaptor {

void Initialize(int numScripts, char* scripts[])
{
	if (Processor == NULL) {
		Processor = vtkCPProcessor::New();
		Processor->Initialize();
	} else {
		Processor->RemoveAllPipelines();
	}

//for(int i=1;i<numScripts;i++)
//  {
	vtkNew<vtkCPPythonScriptPipeline> pipeline;
	pipeline->Initialize(scripts[numScripts - 1]);
	Processor->AddPipeline(pipeline.GetPointer());
//  }
}

void Finalize() {
	if (Processor) {
		Processor->Delete();
		Processor = NULL;
	}
	if (VTKGrid) {
		VTKGrid->Delete();
		VTKGrid = NULL;
	}
}

void CoProcess(
		mesh::Mesh *mesh,
		std::vector<std::vector<eslocal> > &l2g_vec,
		std::vector<std::vector<double> >& prim_solution,
		double time,
		unsigned int timeStep,
		bool lastTimeStep)
{
	vtkNew<vtkCPDataDescription> dataDescription;
	dataDescription->AddInput("input");
	dataDescription->SetTimeData(time, timeStep);
	if (lastTimeStep == true) {
// assume that we want to all the pipelines to execute if it
// is the last time step.
		dataDescription->ForceOutputOn();
	}
	if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
		BuildVTKDataStructures(mesh, l2g_vec, prim_solution);//, decomposition_values);
		dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
		Processor->CoProcess(dataDescription.GetPointer());
	}
}

} // end of Adaptor namespace
