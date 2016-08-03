#include <vtkArrowSource.h>
#include <vtkGlyph2D.h>
#include <vtkGlyph2D.h>
#include <vtkReverseSense.h>
#include <vtkMaskPoints.h>
#include <vtkHedgeHog.h>
#include <vtkBrownianPoints.h>
#include <vtkSphereSource.h>

#include <vtkGenericDataObjectReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDecimatePro.h>
#include <vtkAppendFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkExtractUnstructuredGrid.h>
#include <vtkContourGrid.h>
#include <vtkContourFilter.h>
#include <vtkEnSightWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <string>
#include <vtkPoints.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMPIController.h>

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
#include <vtkMultiProcessController.h>
#include <vtkMPICommunicator.h>
#include <vtkAppendPolyData.h>
#include <vtkMergeCells.h>


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

VTK::VTK(const Mesh &mesh, const std::string &path): Results(mesh, path)
{
	// constructor
}

void VTK::mesh(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shrinkCluster)
{
	std::cout << "SAVE GENERIC VTK DATA\n";
}

void VTK::properties(const Mesh &mesh, const std::string &path, std::vector<Property> properties, double shrinkSubdomain, double shrinkCluster)
{
	std::cout << "SAVE PROPERTIES TO VTK\n";
}

void VTK::fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	std::cout << "SAVE FIX POINTS TO VTK\n";
}

void VTK::corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster)
{
	std::cout << "SAVE CORNERS TO VTK\n";
}

void VTK::store(std::vector<std::vector<double> > &displasment, double shrinkSubdomain, double shrinkCluster)
{
	//TO DO compresion
	std::cout << "Compression: " << config::output::OUTPUT_COMPRESSION << "\n";
//pokus();

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
	vtkNew < vtkPoints > points;
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
			VTKGrid->InsertNextCell(elements[i]->vtkCode(), elements[i]->size(),
					&tmp[0]);
			i++;
		}
		cnt += _coordinates.localSize(part);
	}
//decompositon

	float *decomposition_array = new float[elements.size()];

	counter = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			float part_redefine = part;
			decomposition_array[counter] = part_redefine;
			counter++;
		}
	}

	vtkNew < vtkFloatArray > decomposition;
	decomposition->SetName("decomposition");
	decomposition->SetNumberOfComponents(1);
	decomposition->SetArray(decomposition_array,
			static_cast<vtkIdType>(elements.size()), 1);
	VTKGrid->GetCellData()->AddArray(decomposition.GetPointer());

	float *d_array = new float[elements.size()];
	counter = 0;
	for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
		for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
			float part_redefine = part;
			if (part < (_partPtrs.size() / 2)) {
				d_array[counter] = 1;
			} else {
				d_array[counter] = 0;
			}
			counter++;
		}
	}

	vtkNew < vtkFloatArray > dec;
	dec->SetName("pokus");
	dec->SetNumberOfComponents(1);
	dec->SetArray(d_array, static_cast<vtkIdType>(elements.size()), 0);
	VTKGrid->GetCellData()->AddArray(dec.GetPointer());

//displacement
	int mycounter = 0;
	for (size_t i = 0; i < displasment.size(); i++) {
		mycounter += displasment[i].size();
	}

	double displacement_array[mycounter];

	counter = 0;

	for (size_t i = 0; i < displasment.size(); i++) {
		for (size_t j = 0; j < (displasment[i].size() / 3); j++) {
			displacement_array[3 * counter + 0] = displasment[i][3 * j + 0];
			displacement_array[3 * counter + 1] = displasment[i][3 * j + 1];
			displacement_array[3 * counter + 2] = displasment[i][3 * j + 2];

			counter++;
		}
	}

	vtkNew < vtkDoubleArray > displacement;
	displacement->SetName("displacement");
	displacement->SetNumberOfComponents(3);
	displacement->SetNumberOfTuples(static_cast<vtkIdType>(counter));
	VTKGrid->GetPointData()->AddArray(displacement.GetPointer());

	double* displacementData = displacement_array;
	vtkIdType numTuples = displacement->GetNumberOfTuples();
	for (vtkIdType i = 0, counter = 0; i < numTuples; i++, counter++) {
		double values[3] = { displacementData[i * 3],
				displacementData[i * 3 + 1], displacementData[i * 3 + 2] };
		displacement->SetTypedTuple(counter, values);
	}

//surface and decimation

	std::cout << "Decimation: " << config::output::OUTPUT_DECIMATION << "\n";
	vtkSmartPointer < vtkDecimatePro > decimate = vtkSmartPointer
			< vtkDecimatePro > ::New();
	vtkAppendFilter* ap = vtkAppendFilter::New();
	vtkUnstructuredGrid* decimated = vtkUnstructuredGrid::New();
	vtkGeometryFilter* gf = vtkGeometryFilter::New();
	vtkDataSetSurfaceFilter* dssf = vtkDataSetSurfaceFilter::New();
	vtkTriangleFilter* tf = vtkTriangleFilter::New();

	vtkSmartPointer < vtkContourFilter > cg = vtkSmartPointer < vtkContourFilter
			> ::New();

	if (config::output::OUTPUT_DECIMATION != 0) {
		gf->SetInputData(VTKGrid);
		gf->Update();

		dssf->SetInputData(gf->GetOutput());
		dssf->Update();

		tf->SetInputData(dssf->GetOutput());
		tf->Update();

		decimate->SetInputData(tf->GetOutput());
		decimate->SetTargetReduction(config::output::OUTPUT_DECIMATION);
		decimate->Update();

		ap->AddInputData(decimate->GetOutput());
		ap->Update();

		decimated->ShallowCopy(ap->GetOutput());
	}

//add BlockId
	bool FCD = false;
	if (VTKGrid->GetCellData()) {
		if (VTKGrid->GetCellData()->GetArray("BlockId")) {
			FCD = true;
		}
	}
	if (FCD == false) {
		vtkIntArray *bids = vtkIntArray::New();
		bids->SetName("BlockId");
		for (i = 0; i < VTKGrid->GetNumberOfCells(); i++) {
			bids->InsertNextValue(1);
		}
		VTKGrid->GetCellData()->SetScalars(bids);
	}
	int blockids[2];
	blockids[0] = 1;
	blockids[1] = 0;

//MultiProces controler

	vtkMergeCells* merge = vtkMergeCells::New();
	vtkUnstructuredGrid* ugcase = vtkUnstructuredGrid::New();

	vtkMPIController* controller = vtkMPIController::New();
	controller->Initialize();
//vtkMultiProcessController* controller = vtkMPIController::New();
//controller->GetGlobalController();
	int rank = config::env::MPIrank;
	if (rank != 0) {
		controller->Send(VTKGrid, 0, 1111 + rank);
	}
//controller->Broadcast(h,rank);

//Compresor
	vtkZLibDataCompressor *myZlibCompressor = vtkZLibDataCompressor::New();
	myZlibCompressor->SetCompressionLevel(9);

//vtk
	vtkSmartPointer < vtkGenericDataObjectWriter > writervtk = vtkSmartPointer
			< vtkGenericDataObjectWriter > ::New();
	std::stringstream ss;

//vtu
	vtkSmartPointer < vtkXMLUnstructuredGridWriter > writervtu = vtkSmartPointer
			< vtkXMLUnstructuredGridWriter > ::New();
	std::stringstream sss;

//vtm
	vtkSmartPointer < vtkXMLUnstructuredGridWriter > wvtu = vtkSmartPointer
			< vtkXMLUnstructuredGridWriter > ::New();
	std::stringstream ssss;
	ofstream result;
//ensight
	vtkEnSightWriter *wcase = vtkEnSightWriter::New();

	switch (config::output::OUTPUT_FORMAT) {
	case config::output::OUTPUT_FORMATAlternatives::VTK_LEGACY_FORMAT:
		std::cout << "LEGACY\n";
		ss << _path << ".vtk";
		writervtk->SetFileName(ss.str().c_str());
		if (config::output::OUTPUT_DECIMATION == 0) {
			writervtk->SetInputData(VTKGrid);
		} else {
			writervtk->SetInputData(decimated);
		}
		if (config::output::OUTPUT_COMPRESSION) {
			// writervtk->SetCompressor(myZlibCompressor);
			std::cout << "Compression: " << config::output::OUTPUT_COMPRESSION
					<< "\n";
		}

		writervtk->Write();

		break;

	case config::output::OUTPUT_FORMATAlternatives::VTK_BINARY_FORMAT:
		std::cout << "BINARY\n";
		sss << _path << ".vtu";
		writervtu->SetFileName(sss.str().c_str());
		if (config::output::OUTPUT_DECIMATION == 0) {
			writervtu->SetInputData(VTKGrid);
		} else {
			writervtu->SetInputData(decimated);
		}
		writervtu->SetDataModeToBinary();
		writervtu->Write();
		break;

	case config::output::OUTPUT_FORMATAlternatives::VTK_MULTIBLOCK_FORMAT:
		std::cout << "MULTIBLOCK\n";
		ssss << _path << ".vtu";
		wvtu->SetFileName(ssss.str().c_str());
		if (config::output::OUTPUT_DECIMATION == 0) {
			wvtu->SetInputData(VTKGrid);
		} else {
			wvtu->SetInputData(decimated);
		}
		wvtu->SetDataModeToBinary();
		wvtu->Write();

		MPI_Barrier(MPI_COMM_WORLD);
		int size = config::env::MPIsize;
		if (config::env::MPIrank == 0) {
			result.open("result.vtm");
			result << "<?xml version=\"1.0\"?>\n";
			if (!config::output::OUTPUT_COMPRESSION) {
				result
						<< "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
				std::cout << "Compression: "
						<< config::output::OUTPUT_COMPRESSION << "\n";
			} else {
				result
						<< "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
				std::cout << "Compression: "
						<< config::output::OUTPUT_COMPRESSION << "\n";
			}

			result << " <vtkMultiBlockDataSet>\n";
			for (int i = 0; i < size; i++) {
				result << "  <DataSet index=\"" << i << "\" file=\"result" << i
						<< ".vtu\">\n  </DataSet>\n";
			}
			result << " </vtkMultiBlockDataSet>\n";
			result << "</VTKFile>\n";
			result.close();
		}
		break;

	case config::output::OUTPUT_FORMATAlternatives::ENSIGHT_FORMAT:
		std::cout << "ENSIGHT\n";

		//write ensight format
		wcase->SetFileName("result.case");
		wcase->SetNumberOfBlocks(1);
		wcase->SetBlockIDs(blockids);
		wcase->SetTimeStep(0);
		if (config::env::MPIrank == 0) {
			vtkAppendFilter* app = vtkAppendFilter::New();
			app->AddInputData(VTKGrid);
			for (int i = 1; i < config::env::MPIsize; i++) {
				vtkUnstructuredGrid* h = vtkUnstructuredGrid::New();
				controller->Receive(h, i, 1111 + i);
				app->AddInputData(h);
				std::cout << "prijato\n";
			}
			app->Update();
			ugcase->ShallowCopy(app->GetOutput());
			wcase->SetInputData(ugcase);
			wcase->Write();
			wcase->WriteCaseFile(1);

			//pokus
			vtkSmartPointer < vtkSphereSource > sphereSource = vtkSmartPointer
					< vtkSphereSource > ::New();

			vtkMath::RandomSeed(5070); // for testing
			vtkSmartPointer < vtkBrownianPoints > brownianPoints =
					vtkSmartPointer < vtkBrownianPoints > ::New();
			brownianPoints->SetInputConnection(sphereSource->GetOutputPort());

			vtkSmartPointer < vtkArrowSource > arrowSource = vtkSmartPointer
					< vtkArrowSource > ::New();

			vtkSmartPointer < vtkGlyph3D > glyph3D = vtkSmartPointer
					< vtkGlyph3D > ::New();
			glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
			glyph3D->SetInputConnection(brownianPoints->GetOutputPort());
			glyph3D->SetScaleFactor(.3);

			/*vtkHedgeHog* hh=vtkHedgeHog::New();
			 hh->SetInputData(polydata);
			 hh->SetOutputPointsPrecision(11);
			 hh->Update();*/

			vtkSmartPointer < vtkPolyDataWriter > pw = vtkSmartPointer
					< vtkPolyDataWriter > ::New();
			pw->SetFileName("pokus.vtk");
			pw->SetInputData(glyph3D->GetOutput());
			pw->Write();

			//pokus

		}

		break;
	}
}
