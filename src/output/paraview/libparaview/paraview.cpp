
#include "../paraview.h"

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include <sstream>
#include <vtkGenericDataObjectWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkMPIController.h>
#include <vtkEnSightWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkAppendFilter.h>
#include <vtkGeometryFilter.h>

vtkSmartPointer<vtkCPProcessor> processor;
vtkSmartPointer<vtkUnstructuredGrid> VTKGrid;
vtkNew<vtkCPDataDescription> dataDescription;
double xx = 1;

using namespace espreso::output;

Paraview::Paraview(const Mesh &mesh, const std::string &path): Store(mesh, path)
{
	processor = vtkSmartPointer<vtkCPProcessor>::New();
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
	const unsigned int dofs= displasment[0].size()/ _mesh.coordinates().localSize(0);
	unsigned int size = mycounter;

	mycounter = 0;
	for (size_t i = 0; i < displasment.size(); i++) {
		//note prim_solution.size = number of subdomains. all elements in prim_solution[0]
		for (size_t j = 0; j < (displasment[i].size() / dofs); j++) {
			for(size_t k=0;k<dofs;k++){
				displacement_array[dofs*mycounter+k]=displasment[i][dofs*j+k];

			}
			mycounter++;

		}
	}

	if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0) {
		// displacement array
		vtkNew<vtkDoubleArray> displacement;
		displacement->SetName("displacement");
		displacement->SetNumberOfComponents(dofs);
		displacement->SetNumberOfTuples(static_cast<vtkIdType>(mycounter));
		VTKGrid->GetPointData()->AddArray(displacement.GetPointer());
	}
	vtkDoubleArray* displacement = vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray("displacement"));
	double* displacementData = displacement_array;
	vtkIdType numTuples = displacement->GetNumberOfTuples();
	for (vtkIdType i = 0, counter = 0; i < numTuples; i++, counter++) {
		double *values=new double[dofs];
		for(size_t k=0;k<dofs;k++){
			values[k]=displacementData[i*dofs+k];
		}
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
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<eslocal> &_partPtrs = _mesh.getPartition();

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
	decomposition->SetArray(decomposition_array, static_cast<vtkIdType>(elements.size()), numb);
	numb++;
	VTKGrid->GetCellData()->AddArray(decomposition.GetPointer());


}

void Paraview::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<eslocal> &partition = _mesh.getPartition();
	std::vector<std::vector<double> > selection(_mesh.parts());
	std::vector<std::vector<double> > values(_mesh.parts());

	vtkSmartPointer<vtkUnstructuredGrid> prof = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> po = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> sv = vtkSmartPointer<vtkFloatArray>::New();
	sv->SetName(name.c_str());
	sv->SetNumberOfComponents(3);

	vtkIdType it = 0;
	vtkSmartPointer<vtkCellArray> ver = vtkSmartPointer<vtkCellArray>::New();

	if (eType == ElementType::NODES) {
		int tupl=0;
		for (size_t n = 0; n < _mesh.nodes().size(); n++) {
			const Element *node = _mesh.nodes()[n];
			for (size_t p = 0; p < properties.size(); p++) {
				for (size_t d = 0; d < node->domains().size(); d++) {
					size_t domain = node->domains()[d];
					selection[domain].push_back(node->settings().isSet(properties[p]) ? 1 : 0);
					values[domain].push_back(node->settings(properties[p]).back()->evaluate(n));
					tupl++;
				}
			}
		}

		vtkNew < vtkDoubleArray > selections;
		std::string names=name+"IsSet";
		selections->SetName(names.c_str());
		selections->SetNumberOfComponents(properties.size());
		selections->SetNumberOfTuples(static_cast<vtkIdType>(tupl));
		VTKGrid->GetPointData()->AddArray(selections.GetPointer());

		int c=0;
		for(size_t p=0;p<selection.size();p++){
			for (size_t i = 0; i < selection[p].size()/properties.size(); i++) {
				double* select=new double[properties.size()];
				//std::vector<double> select(properties.size());
				for (size_t k = 0; k < properties.size(); k++) {
					select[k] = {selection[p][i*properties.size()+k]};
				}
				selections->SetTypedTuple(c, select);
				c++;
				delete [] select;
			}
		}

		vtkNew < vtkDoubleArray > valuess;
		std::string namess=name+"FixedValue";
		valuess->SetName(namess.c_str());
		valuess->SetNumberOfComponents(properties.size());
		valuess->SetNumberOfTuples(static_cast<vtkIdType>(tupl));
		VTKGrid->GetPointData()->AddArray(valuess.GetPointer());

		c=0;
		for(size_t p=0;p<values.size();p++){
			for (size_t i = 0; i < values[p].size()/properties.size(); i++) {
				double* select=new double[properties.size()];
				for (size_t k = 0; k < properties.size(); k++) {
					select[k] = {values[p][i*properties.size()+k]};
				}
				valuess->SetTypedTuple(c, select);
				c++;
				delete [] select;
			}
		}
	} else if (eType == ElementType::ELEMENTS) {
		float* hvalues=new float[_mesh.elements().size()*properties.size()];

		for (size_t e = 0; e < _mesh.elements().size(); e++) {
			const Element *element = _mesh.elements()[e];
			for (size_t p = 0; p < properties.size(); p++) {
				size_t domain = element->domains()[0];
				double value = 0;
				for (size_t n = 0; n < element->nodes(); n++) {
					value += element->settings(properties[p]).back()->evaluate(element->node(n));
				}
				values[domain].push_back(value / element->nodes());
			}
		}
		values.push_back(std::vector<double>(_mesh.faces().size() + _mesh.edges().size()));

		int ite=0;
		for(size_t i=0;i<values.size();i++){
			for(size_t j=0;j<values[i].size()/properties.size();j++){
				for(size_t d=0;d<properties.size();d++){
					hvalues[ite]=values[i][properties.size()*j+d];
					ite++;
				}
			}
		}
		vtkNew < vtkFloatArray > value;
		value->SetName(name.c_str());
		value->SetNumberOfComponents(properties.size());
		value->SetArray(hvalues,static_cast<vtkIdType>(ite), numb);
		numb++;
		VTKGrid->GetCellData()->AddArray(value.GetPointer());
	}

}

void Paraview::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<eslocal> &_partPtrs = _mesh.getPartition();

	//Values
	int counter;
	if (eType == ElementType::NODES) {
		int mycounter = 0;
		for (size_t i = 0; i < values.size(); i++) {
			mycounter += values[i].size();
		}

		const unsigned int dofs = values[0].size() / _mesh.coordinates().localSize(0);
		double value_array[mycounter];
		counter = 0;

		for (size_t i = 0; i < values.size(); i++) {
			for (size_t j = 0; j < (values[i].size() / dofs); j++) {
				for (int k = 0; k < dofs; k++) {
					value_array[dofs * counter + k] = values[i][dofs * j + k];
				}
				counter++;
			}
		}

		vtkNew<vtkDoubleArray> value;
		value->SetName(name.c_str());
		value->SetNumberOfComponents(dofs);
		value->SetNumberOfTuples(static_cast<vtkIdType>(counter));
		VTKGrid->GetPointData()->AddArray(value.GetPointer());

		double* valueData = value_array;
		vtkIdType numTuples = value->GetNumberOfTuples();
		for (vtkIdType i = 0, counter = 0; i < numTuples; i++, counter++) {
			double values[dofs];
			for (int k = 0; k < dofs; k++) {
				values[k] = {valueData[i*dofs+k]};
			}
			value->SetTypedTuple(counter, values);
		}
	} else if (eType == ElementType::ELEMENTS) {
		float *value_array = new float[elements.size()];

		counter = 0;
		for (size_t part = 0; part + 1 < _partPtrs.size(); part++) {
			for (eslocal i = 0; i < _partPtrs[part + 1] - _partPtrs[part]; i++) {
				float part_redefine = part;
				value_array[counter] = part_redefine;
				counter++;
			}
		}

		vtkNew<vtkFloatArray> value;
		value->SetName(name.c_str());
		value->SetNumberOfComponents(1);
		value->SetArray(value_array, static_cast<vtkIdType>(elements.size()), numb);
		numb++;
		VTKGrid->GetCellData()->AddArray(value.GetPointer());
	}
}

void Paraview::finalize()
{
	switch (config::output::OUTPUT_FORMAT) {
		case config::output::OUTPUT_FORMATAlternatives::VTK_LEGACY_FORMAT:{
			vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
			std::stringstream ss;
			ss<<"result"<<config::env::MPIrank<<".vtk";
			writervtk->SetFileName(ss.str().c_str());
			writervtk->SetInputData(VTKGrid);
			writervtk->Write();
		}
		break;
		case config::output::OUTPUT_FORMATAlternatives::VTK_BINARY_FORMAT:{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			std::stringstream ss;
			ss<<"result"<<config::env::MPIrank<<".vtu";
			writervtu->SetFileName(ss.str().c_str());
			writervtu->SetInputData(VTKGrid);
			writervtu->SetDataModeToBinary();
			writervtu->Write();
		}
		break;
		case config::output::OUTPUT_FORMATAlternatives::VTK_MULTIBLOCK_FORMAT:{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			std::stringstream ss;
			ss<<"result"<<config::env::MPIrank<<".vtu";
			wvtu->SetFileName(ss.str().c_str());
			wvtu->SetInputData(VTKGrid);
			wvtu->SetDataModeToBinary();
			wvtu->Write();
		}
		break;
		case config::output::OUTPUT_FORMATAlternatives::ENSIGHT_FORMAT: {
			vtkMPIController* controller=vtkMPIController::New();
			vtkSmartPointer<vtkEnSightWriter> wcase = vtkSmartPointer<vtkEnSightWriter>::New();
			vtkSmartPointer<vtkUnstructuredGrid> ugcase = vtkSmartPointer<vtkUnstructuredGrid>::New();
			bool FCD = false;
			if (VTKGrid->GetCellData()) {
				if (VTKGrid->GetCellData()->GetArray("BlockId")) {
					FCD = true;
				}
			}
			if (FCD == false) {
				vtkSmartPointer<vtkIntArray> bids = vtkSmartPointer<vtkIntArray>::New();
				bids->SetName("BlockId");
				for (int i = 0; i < VTKGrid->GetNumberOfCells(); i++) {
					bids->InsertNextValue(1);
				}
				VTKGrid->GetCellData()->SetScalars(bids);
			}
			int blockids[2];
			blockids[0] = 1;
			blockids[1] = 0;

			if (config::env::MPIrank != 0) {
				controller->Send(VTKGrid, 0, 1111 + config::env::MPIrank);
			}

			if (config::env::MPIrank == 0) {
				wcase->SetFileName("result");
				wcase->SetNumberOfBlocks(1);
				wcase->SetBlockIDs(blockids);
				wcase->SetTimeStep(0);
				vtkSmartPointer<vtkAppendFilter> app = vtkSmartPointer<vtkAppendFilter>::New();
				app->AddInputData(VTKGrid);
				for (int i = 1; i < config::env::MPIsize; i++) {
					vtkSmartPointer<vtkUnstructuredGrid> h = vtkSmartPointer<vtkUnstructuredGrid>::New();
					controller->Receive(h, i, 1111 + i);
					app->AddInputData(h);
				}
				app->Update();
				vtkSmartPointer<vtkGeometryFilter> gf = vtkSmartPointer<vtkGeometryFilter>::New();
				gf->SetInputData(app->GetOutput());
				gf->Update();
				vtkSmartPointer<vtkCleanPolyData> cpd = vtkSmartPointer<vtkCleanPolyData>::New();
				cpd->SetInputData(gf->GetOutput());
				cpd->Update();
				vtkSmartPointer<vtkAppendFilter> apc = vtkSmartPointer<vtkAppendFilter>::New();
				apc->SetInputData(cpd->GetOutput());
				apc->Update();
				ugcase->ShallowCopy(apc->GetOutput());
				wcase->SetInputData(ugcase);
				wcase->Write();
				wcase->WriteCaseFile(1);
			}
			controller->Delete();
		}
		break;
		}
}
