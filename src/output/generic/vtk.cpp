#include <vtkGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <string>
#include <vtkPoints.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkObjectFactory.h>
#include <vtkZLibDataCompressor.h>
#include <vtk_zlib.h>

#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <iostream>
#include <fstream>

#include "esmesh.h"
#include "esoutput.h"

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

using namespace espreso::output;

vtkUnstructuredGrid* VTKGrid;

Generic::Generic(const Mesh &mesh, const std::string &path): ResultStore(mesh, path)
{
	// constructor
}

void Generic::store(double shrinkSubdomain, double shrinkCluster)
{
  std::cout << "SAVE GENERIC VTK DATA\n";
}

void Generic::store(std::vector<std::vector<double> > &displasment, double shrinkSubdomain, double shrinkCluster)
{

	const std::vector<Element*> &elements = _mesh.getElements();
	const std::vector<eslocal> &_partPtrs = _mesh.getPartition();
	
	VTKGrid = vtkUnstructuredGrid::New();        

	size_t n_nodsClust = 0;
	for (size_t iEl = 0; iEl < elements.size(); iEl++) {
		n_nodsClust += elements[iEl]->size();
	}
	size_t cnt = 0, n_points = 0;
	for (size_t d = 0; d < _mesh.parts(); d++) {
		n_points += _mesh.coordinates().localSize(d);
	}

	double shrinking = 0.90;
	const Coordinates &_coordinates = _mesh.coordinates();
	double *coord_xyz = new double[n_points * 3];

	int counter = 0;
	for (size_t d = 0; d < _mesh.parts(); d++) {
		Point center;
		for (size_t c = 0; c < _coordinates.localSize(d); c++) {
			center += _coordinates.get(c, d);
		}
		center /= _coordinates.localSize(d);

		for (size_t i = 0; i < _coordinates.localSize(d); i++) {
			Point xyz = _coordinates.get(i, d);
			xyz = center + (xyz - center) * shrinking;
			coord_xyz[3 * counter + 0] = xyz.x;
			coord_xyz[3 * counter + 1] = xyz.y;
			coord_xyz[3 * counter + 2] = xyz.z;
			counter++;
		}
	}

	vtkNew < vtkDoubleArray > pointArray;
	pointArray->SetNumberOfComponents(3);

	size_t numpoints = n_points * 3;
	pointArray->SetArray(coord_xyz, numpoints, 1);
	vtkNew <vtkPoints > points;
	points->SetData(pointArray.GetPointer());
	VTKGrid->SetPoints(points.GetPointer());

	VTKGrid->Allocate(static_cast<vtkIdType>(n_nodsClust));
	vtkIdType tmp[100]; //max number of  node

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
		cnt += _coordinates.localSize(part);
	}
	//decompositon	
	
	float *decomposition_array = new float[elements.size()];
	
	counter = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			float part_redefine =  part;
			decomposition_array[counter] = part_redefine;
			counter++;
		}
	}
	
	vtkNew<vtkFloatArray> decomposition;
	decomposition->SetName("decomposition");
	decomposition->SetNumberOfComponents(1);	

	//vtkFloatArray* decomposition = vtkFloatArray::SafeDownCast(VTKGrid->GetCellData()->GetArray("decomposition"));
	decomposition->SetArray(decomposition_array, static_cast<vtkIdType>(elements.size()), 1);
	VTKGrid->GetCellData()->AddArray(decomposition.GetPointer());

	//displacement
	int mycounter=0;
	for (size_t i = 0; i < displasment.size(); i++) {
		mycounter += displasment[i].size();
	}
	
	double displacement_array[mycounter];

	counter=0;

	for (size_t i = 0; i < displasment.size(); i++) {
	  for (size_t j = 0; j < (displasment[i].size() / 3); j++) {
			displacement_array[3 * counter + 0] = displasment[i][3 * j + 0];
			displacement_array[3 * counter + 1] = displasment[i][3 * j + 1];
			displacement_array[3 * counter + 2] = displasment[i][3 * j + 2];

			counter++;
		}
	}
	
	vtkNew<vtkDoubleArray> displacement;
	displacement->SetName("displacement");
	displacement->SetNumberOfComponents(3);
	displacement->SetNumberOfTuples(static_cast<vtkIdType>(counter));
	VTKGrid->GetPointData()->AddArray(displacement.GetPointer());

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


	//writer vtu
	

	vtkZLibDataCompressor *myZlibCompressor=vtkZLibDataCompressor::New();
	myZlibCompressor->SetCompressionLevel(25);

	vtkSmartPointer< vtkXMLUnstructuredGridWriter > writer = vtkSmartPointer< vtkXMLUnstructuredGridWriter > ::New();
	std::stringstream ss;
	ss << _path << ".vtu";
	writer->SetFileName(ss.str().c_str());
	writer->SetInputData(VTKGrid);
	writer->SetDataModeToBinary();
	writer->Write();

        
	//mbds->SetBlock(rank,writer->GetInput());

	config::env::MPIrank;
	config::env::MPIsize;
	//writer vtm
	MPI_Barrier(MPI_COMM_WORLD);
	int size=config::env::MPIsize;
	ofstream result;
	result.open("result.vtm");
	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	result<<" <vtkMultiBlockDataSet>\n";
	for(int i=0 ; i < size ; i++){
	  result<<"  <DataSet index=\""<< i <<"\" file=\"result"<< i <<".vtu\">\n  </DataSet>\n";
	}
	result<<" </vtkMultiBlockDataSet>\n";
	result<<"</VTKFile>\n";
	result.close();

	/*vtkMultiBlockDataSet *mbds=vtkMultiBlockDataSet::New();
	vtkXMLMultiBlockDataWriter *vtmwriter=vtkXMLMultiBlockDataWriter::New();	
	vtmwriter->SetFileName("res.vtm");
	
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0;i<config::env::MPIsize;i++){
	  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	  std::stringstream sss;
	  sss << "result"<<i << ".vtu";
	  reader->SetFileName(sss.str().c_str());
	  reader->Update();

	  mbds->SetBlock(i,reader->GetOutput());
	}
	vtmwriter->SetInputData(mbds);	  
	vtmwriter->SetCompressor(myZlibCompressor);
	vtmwriter->Write();*/
	
	/*for(int i=0;i<config::env::MPIsize;i++)
	  {
	    if(config::env::MPIrank==0){
	      mbds->SetBlock(0,VTKGrid);
	      vtmwriter->Update();
	    }
            else if(config::env::MPIrank==i){
	      mbds->SetBlock(i,VTKGrid);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
	  vtmwriter->SetInputData(mbds);
	  vtmwriter->SetCompressor(myZlibCompressor);
	  vtmwriter->Write();
	}
	if(rank==1){
	  vtmwriter->SetInputData(mbds);
	  vtmwriter->SetCompressor(myZlibCompressor);
	  vtmwriter->Update();
	  }*/
	
	std::cout << "SAVE GENERIC VTK RESULT XX\n";
}











/*int main(int argc, char *argv[])
{
  std::string inputFilename=argv[1];
  std::string outputFilename=argv[2];

  //input
  vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();

  //output
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(outputFilename.c_str());
  writer->SetInputConnection(reader->GetOutputPort());
  writer->Update();

  //mapper a actor
  vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInputConnection(reader->GetOutputPort());

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  //renderer
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(1,1,1);

  renderWindow->Render();
  renderWindowInteractor->Start();


  return EXIT_SUCCESS;
}*/
