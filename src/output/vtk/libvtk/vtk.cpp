#include <vtkArrowSource.h>
#include <vtkGlyph2D.h>
#include <vtkGlyph2D.h>
#include <vtkReverseSense.h>
#include <vtkMaskPoints.h>
#include <vtkHedgeHog.h>
#include <vtkBrownianPoints.h>
#include <vtkSphereSource.h>
#include <vtkStructuredGrid.h>
#include <vtkLineSource.h>
#include <vtkLine.h>
#include <vtkTriangle.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkXMLCompositeDataWriter.h>
#include <vtkCleanPolyData.h>

#include <vtkGenericDataObjectReader.h>
#include <vtkEnSightReader.h>
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
#include <vtkAppendPolyData.h>
#include <vtkMergeCells.h>

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

#include <vtkMultiProcessController.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>

#include "../vtk.h"

using namespace espreso::store;

VTK::VTK(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: Store(output, mesh, path), _lastData(ElementType::ELEMENTS)
{
	computeCenters();
}

void VTK::storeGeometry(size_t timeStep)
{
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<eslocal> &_partPtrs = _mesh.getPartition();

	vtkSmartPointer<vtkUnstructuredGrid> VTKGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

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
			//xyz=_sCenters[d] + (xyz - _sCenters[d]) * _shrinkSubdomain;
			xyz=shrink(xyz,d);
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


	stringstream ss;


	//Writers
	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT: {
		vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
		ss << _path << environment->MPIrank << ".vtk";
		writervtk->SetFileName(ss.str().c_str());
		writervtk->SetInputData(VTKGrid);
		writervtk->Write();
	}
		break;
	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		ss << _path << environment->MPIrank << ".vtu";
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(VTKGrid);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
		break;
	}
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		std::ofstream result;

		ss << _path << environment->MPIrank << ".vtu";
		wvtu->SetFileName(ss.str().c_str());
		wvtu->SetInputData(VTKGrid);
		wvtu->SetDataModeToBinary();
		wvtu->Write();

		int size = environment->MPIsize;
		if (environment->MPIrank == 0) {
			stringstream sss;
			sss << _path << ".vtm";
			result.open(sss.str().c_str());
			result << "<?xml version=\"1.0\"?>\n";
			if (!_output.compression) {
				result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";

			} else {
				result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";

			}

			result << " <vtkMultiBlockDataSet>\n";
			for (int i = 0; i < size; i++) {
				result << "  <DataSet index=\"" << i << "\" file=\"" << _path << i << ".vtu\">\n  </DataSet>\n";
			}
			result << " </vtkMultiBlockDataSet>\n";
			result << "</VTKFile>\n";
			result.close();
		}
	}
	break;
	case OUTPUT_FORMAT::ENSIGHT_FORMAT: {
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(VTKGrid);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
	}
	break;
	}
	delete[] decomposition_array;
}

void VTK::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<eslocal> &partition = _mesh.getPartition();
	std::vector<std::vector<double> > selection(_mesh.parts());
	std::vector<std::vector<double> > values(_mesh.parts());

	vtkSmartPointer<vtkPolyData> pro = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkUnstructuredGrid> prof = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> po = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> sv = vtkSmartPointer<vtkFloatArray>::New();
	sv->SetName(name.c_str());
	sv->SetNumberOfComponents(3);

	stringstream ss;


	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
		vtkSmartPointer<vtkGenericDataObjectReader> rvtk = vtkSmartPointer<vtkGenericDataObjectReader>::New();
		ss<<_path<<environment->MPIrank<<".vtk";
		rvtk->SetFileName(ss.str().c_str());
		rvtk->Update();
		prof->ShallowCopy(rvtk->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rvtu = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rvtu->SetFileName(ss.str().c_str());
		rvtu->Update();
		prof->ShallowCopy(rvtu->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rvtu = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rvtu->SetFileName(ss.str().c_str());
		rvtu->Update();
		prof->ShallowCopy(rvtu->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::ENSIGHT_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rcase = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rcase->SetFileName(ss.str().c_str());
		rcase->Update();
		prof->ShallowCopy(rcase->GetOutput());
	}
	break;
	}

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
		string names=name+"IsSet";
		selections->SetName(names.c_str());
		selections->SetNumberOfComponents(properties.size());
		selections->SetNumberOfTuples(static_cast<vtkIdType>(tupl));
		prof->GetPointData()->AddArray(selections.GetPointer());

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
		string namess=name+"FixedValue";
		valuess->SetName(namess.c_str());
		valuess->SetNumberOfComponents(properties.size());
		valuess->SetNumberOfTuples(static_cast<vtkIdType>(tupl));
		prof->GetPointData()->AddArray(valuess.GetPointer());

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
		prof->GetCellData()->AddArray(value.GetPointer());
		delete[] hvalues;
	}

	stringstream sss;
	//Writers
	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
		vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
		writervtk->SetFileName(ss.str().c_str());
		writervtk->SetInputData(prof);
		writervtk->Write();
	}
	break;
	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(prof);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
	}
	break;
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		std::ofstream result;
		wvtu->SetFileName(ss.str().c_str());
		wvtu->SetInputData(prof);
		wvtu->SetDataModeToBinary();
		wvtu->Write();
	}
	break;
	case OUTPUT_FORMAT::ENSIGHT_FORMAT: {
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(prof);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
		}
		break;
	}
}

void VTK::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<eslocal> &_partPtrs = _mesh.getPartition();

	vtkSmartPointer<vtkUnstructuredGrid> VTKGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	stringstream ss;
	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
		vtkSmartPointer<vtkGenericDataObjectReader> rvtk = vtkSmartPointer<vtkGenericDataObjectReader>::New();
		ss<<_path<<environment->MPIrank<<".vtk";
		rvtk->SetFileName(ss.str().c_str());
		rvtk->Update();
		VTKGrid->ShallowCopy(rvtk->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rvtu = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rvtu->SetFileName(ss.str().c_str());
		rvtu->Update();
		VTKGrid->ShallowCopy(rvtu->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rvtm = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rvtm->SetFileName(ss.str().c_str());
		rvtm->Update();
		VTKGrid->ShallowCopy(rvtm->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::ENSIGHT_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rcase = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rcase->SetFileName(ss.str().c_str());
		rcase->Update();
		VTKGrid->ShallowCopy(rcase->GetOutput());
	}
	break;
	}

	//Values
	int counter;
	if (eType == ElementType::NODES) {
		int mycounter = 0;
		for (size_t i = 0; i < values.size(); i++) {
			mycounter += values[i].size();
		}

		const unsigned int dofs = values[0].size() / _mesh.coordinates().localSize(0);
		double *value_array=new double[mycounter];
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
		//delete[] value_array;
		delete[] valueData;
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
		delete[] value_array;
	}


	//Writers
	stringstream sss;
	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
		vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
		writervtk->SetFileName(ss.str().c_str());
		writervtk->SetInputData(VTKGrid);
		writervtk->Write();
	}
	break;

	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(VTKGrid);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
	}
	break;

	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		wvtu->SetFileName(ss.str().c_str());
		wvtu->SetInputData(VTKGrid);
		wvtu->SetDataModeToBinary();
		wvtu->Write();
	}
	break;

	case OUTPUT_FORMAT::ENSIGHT_FORMAT: {
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(VTKGrid);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
	}
	break;
	}
}

void VTK::mesh(const OutputConfiguration &output, const Mesh &mesh, const std::string &path, ElementType eType)
{
	std::cout<<"Ukladam mesh..."<<path<<std::endl;
	const std::vector<Element*> &elements = mesh.elements();
	const std::vector<eslocal> &_partPtrs = mesh.getPartition();
	VTK* help = new VTK(output, mesh, path);
	vtkSmartPointer<vtkUnstructuredGrid> VTKGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	size_t n_nodsClust = 0;
	for (size_t iEl = 0; iEl < elements.size(); iEl++) {
		n_nodsClust += elements[iEl]->nodes();
	}
	size_t cnt = 0, n_points = 0;
	for (size_t d = 0; d < mesh.parts(); d++) {
		n_points += mesh.coordinates().localSize(d);
	}

	//Points
	int counter = 0;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t d = 0; d < mesh.parts(); d++) {
		for (size_t i = 0; i < mesh.coordinates().localSize(d); i++) {
			Point xyz = mesh.coordinates().get(i, d);
			//xyz=_sCenters[d] + (xyz - _sCenters[d]) * _shrinkSubdomain;
			xyz=help->shrink(xyz,d);
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
			for (size_t j = 0; j < mesh.elements()[i]->nodes(); j++) {
				tmp[j] = mesh.coordinates().localIndex(elements[i]->node(j), part) + cnt;
			}
			int code = mesh.elements()[i]->vtkCode();
			VTKGrid->InsertNextCell(code, elements[i]->nodes(), &tmp[0]);
			i++;
		}
		cnt += mesh.coordinates().localSize(part);
	}

	stringstream ss;

	//Writers
	switch (output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
		vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
		ss << path << environment->MPIrank << ".vtk";
		writervtk->SetFileName(ss.str().c_str());
		writervtk->SetInputData(VTKGrid);
		writervtk->Write();
	}
	break;
	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		ss << path << environment->MPIrank << ".vtu";
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(VTKGrid);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
		break;
	}
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		std::ofstream result;

		ss << path << environment->MPIrank << ".vtu";
		wvtu->SetFileName(ss.str().c_str());
		wvtu->SetInputData(VTKGrid);
		wvtu->SetDataModeToBinary();
		wvtu->Write();

		int size = environment->MPIsize;
		if (environment->MPIrank == 0) {
			stringstream sss;
			sss << path << ".vtm";
			result.open(sss.str().c_str());
			result << "<?xml version=\"1.0\"?>\n";
			if (!output.compression) {
				result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";

			} else {
				result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";

			}

			result << " <vtkMultiBlockDataSet>\n";
			for (int i = 0; i < size; i++) {
				result << "  <DataSet index=\"" << i << "\" file=\"" << path << i << ".vtu\">\n  </DataSet>\n";
			}
			result << " </vtkMultiBlockDataSet>\n";
			result << "</VTKFile>\n";
			result.close();
		}
	}
	break;
	case OUTPUT_FORMAT::ENSIGHT_FORMAT: {
		ss << path;
		vtkMPIController* controller=vtkMPIController::New();
		controller->Initialize();
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

		if (environment->MPIrank != 0) {
			controller->Send(VTKGrid, 0, 1111 + environment->MPIrank);
		}

		if (environment->MPIrank == 0) {
			wcase->SetFileName(ss.str().c_str());
			wcase->SetNumberOfBlocks(1);
			wcase->SetBlockIDs(blockids);
			wcase->SetTimeStep(0);
			vtkSmartPointer<vtkAppendFilter> app = vtkSmartPointer<vtkAppendFilter>::New();
			app->AddInputData(VTKGrid);
			for (int i = 1; i < environment->MPIsize; i++) {
				vtkSmartPointer<vtkUnstructuredGrid> h = vtkSmartPointer<vtkUnstructuredGrid>::New();
				controller->Receive(h, i, 1111 + i);
				app->AddInputData(h);
			}
			app->Update();
			//std::cout<<app->GetOutput()->GetNumberOfPoints()<<std::endl;
			vtkSmartPointer<vtkGeometryFilter> gf = vtkSmartPointer<vtkGeometryFilter>::New();
			gf->SetInputData(app->GetOutput());
			gf->Update();
			vtkSmartPointer<vtkCleanPolyData> cpd = vtkSmartPointer<vtkCleanPolyData>::New();
			cpd->SetInputData(gf->GetOutput());
			cpd->Update();
			vtkSmartPointer<vtkAppendFilter> apc = vtkSmartPointer<vtkAppendFilter>::New();
			apc->SetInputData(cpd->GetOutput());
			apc->Update();
			//std::cout<<apc->GetOutput()->GetNumberOfPoints()<<std::endl;
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

void VTK::fixPoints(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
{
	vtkSmartPointer<vtkUnstructuredGrid> fp = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//std::vector<std::vector<eslocal> > fixPoints(mesh.parts());
	const Coordinates &_coordinates = mesh.coordinates();
	const std::vector<eslocal> &_partPtrs = mesh.getPartition();
	std::vector<Element*> fixPoints;
	size_t parts = _coordinates.parts();
	double shrinking = 0.90;

	//points
	/*for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t i = 0; i < mesh.fixPoints(p).size(); i++) {
			fixPoints[p].push_back(mesh.fixPoints(p)[i]->node(0));
		}
	}*/
	for (size_t p = 0; p < mesh.parts(); p++) {
		fixPoints.insert(fixPoints.end(), mesh.fixPoints(p).begin(), mesh.fixPoints(p).end());
	}
	std::sort(fixPoints.begin(),fixPoints.end());
	Esutils::removeDuplicity(fixPoints);

	int counter = 0;
	vtkSmartPointer<vtkPoints> point=vtkSmartPointer<vtkPoints>::New();
	for (size_t d = 0; d < fixPoints.size(); d++) {
		for (size_t i = 0; i < fixPoints[d]->domains().size(); i++) {
			Point xyz = _coordinates[fixPoints[d]->node(0)];

			point->InsertNextPoint(xyz.x,xyz.y,xyz.z);
			counter++;
		}
	}
	fp->SetPoints(point);
	//cells
	size_t offset = 0;
	size_t cnt = 0;
	size_t i = 0;
	vtkIdType tmp[100]; //max number of  nodevtkIdType

	for (size_t p = 0; p < fixPoints.size(); p++) {
		for (size_t j = 0; j < fixPoints[i]->domains().size(); j++) {
			tmp[j] = cnt + j;
		}
		fp->InsertNextCell(2, fixPoints[p]->domains().size(), &tmp[0]);

		i++;
		cnt += fixPoints[p]->domains().size();
	}
	//decomposition
	float *decomposition_array = new float[fixPoints.size()];

	for (size_t p = 0; p < fixPoints.size(); p++) {
		decomposition_array[p] = p;

	}

	vtkNew<vtkFloatArray> decomposition;
	decomposition->SetName("decomposition");
	decomposition->SetNumberOfComponents(1);
	decomposition->SetArray(decomposition_array, static_cast<vtkIdType>(fixPoints.size()), 1);
	fp->GetCellData()->AddArray(decomposition.GetPointer());


	//vtk
	vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	std::stringstream ss;

	//vtu
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	std::stringstream sss;

	//vtm
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	std::stringstream ssss;
	ofstream result;
	//ensight
	vtkSmartPointer<vtkEnSightWriter> wcase = vtkSmartPointer<vtkEnSightWriter>::New();

	switch (output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:
		ss << path << environment->MPIrank << ".vtk";
		writervtk->SetFileName(ss.str().c_str());
		writervtk->SetInputData(fp);
		writervtk->Write();
		break;

	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:
		sss << path << environment->MPIrank << ".vtu";
		writervtu->SetFileName(sss.str().c_str());
		writervtu->SetInputData(fp);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
		break;

	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:
		ssss << path << environment->MPIrank << ".vtu";
		wvtu->SetFileName(ssss.str().c_str());
		wvtu->SetInputData(fp);
		wvtu->SetDataModeToBinary();
		wvtu->Write();

		MPI_Barrier(MPI_COMM_WORLD);
		int size = environment->MPIsize;
		if (environment->MPIrank == 0) {
			result.open("meshFP_result.vtm");
			result << "<?xml version=\"1.0\"?>\n";
			result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
			result << " <vtkMultiBlockDataSet>\n";
			for (int i = 0; i < size; i++) {
				result << "  <DataSet index=\"" << i
						<< "\" file=\"meshFixPoints" << i
						<< ".vtu\">\n  </DataSet>\n";
			}
			result << " </vtkMultiBlockDataSet>\n";
			result << "</VTKFile>\n";
			result.close();
		}
		break;

	case OUTPUT_FORMAT::ENSIGHT_FORMAT:
		vtkSmartPointer<vtkUnstructuredGrid> ugcase = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkMPIController* controller = vtkMPIController::New();
		controller->Initialize();
		//add BlockId
		bool FCD = false;
		if (fp->GetCellData()) {
			if (fp->GetCellData()->GetArray("BlockId")) {
				FCD = true;
			}
		}
		if (FCD == false) {
			vtkSmartPointer<vtkIntArray> bids = vtkSmartPointer<vtkIntArray>::New();
			bids->SetName("BlockId");
			for (i = 0; i < fp->GetNumberOfCells(); i++) {
				bids->InsertNextValue(1);
			}
			fp->GetCellData()->SetScalars(bids);
		}
			int blockids[2];
			blockids[0] = 1;
			blockids[1] = 0;
		//write ensight format
		wcase->SetFileName("meshFP_result.case");
		wcase->SetNumberOfBlocks(1);
		wcase->SetBlockIDs(blockids);
		wcase->SetTimeStep(0);
		if(environment->MPIrank!=0){
			int rank = environment->MPIrank;
			controller->Send(fp, 0, 1111 + rank);
		}
		if (environment->MPIrank == 0) {
			vtkSmartPointer<vtkAppendFilter> app = vtkSmartPointer<vtkAppendFilter>::New();
			app->AddInputData(fp);
			for (int i = 1; i < environment->MPIsize; i++) {
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
		break;
	}
	delete[] decomposition_array;
}

void VTK::corners(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
{
	vtkSmartPointer<vtkUnstructuredGrid> c = vtkSmartPointer<vtkUnstructuredGrid>::New();
	const Coordinates &_coordinates = mesh.coordinates();
	const std::vector<eslocal> &_partPtrs = mesh.getPartition();
	const std::vector<Element*> &elements = mesh.elements();
	size_t parts = _coordinates.parts();
	std::vector<std::vector<eslocal> > corners(mesh.parts());

	for (size_t i = 0; i < mesh.corners().size(); i++) {
		for (size_t d = 0; d < mesh.corners()[i]->domains().size(); d++) {
			size_t p = mesh.corners()[i]->domains()[d];
			corners[p].push_back(mesh.corners()[i]->node(0));
		}
	}

	size_t n_nodsClust = 0;
	for (size_t iEl = 0; iEl < elements.size(); iEl++) {
		n_nodsClust += elements[iEl]->nodes();
	}
	vtkSmartPointer<vtkPoints> point=vtkSmartPointer<vtkPoints>::New();
	int counter = 0;
	for (size_t d = 0; d < mesh.parts(); d++) {
		for (size_t i = 0; i < corners[d].size(); i++) {
			Point xyz = mesh.coordinates()[corners[d][i]];
			point->InsertNextPoint(xyz.x,xyz.y, xyz.z);
			counter++;
		}
	}
	c->SetPoints(point);

	//cells
	size_t offset = 0;
	size_t cnt = 0;
	size_t i = 0;
	vtkIdType tmp[100]; //max number of  nodevtkIdType

	for (size_t p = 0; p < corners.size(); p++) {
		for (size_t j = 0; j < corners[i].size(); j++) {
			tmp[j] = cnt + j;
		}
		c->InsertNextCell(2, corners[p].size(), &tmp[0]);

		i++;
		cnt += corners[p].size();
	}
	//decomposition
	float *decomposition_array = new float[corners.size()];

	for (size_t p = 0; p < corners.size(); p++) {
		decomposition_array[p] = p;

	}

	vtkNew<vtkFloatArray> decomposition;
	decomposition->SetName("decomposition");
	decomposition->SetNumberOfComponents(1);
	decomposition->SetArray(decomposition_array, static_cast<vtkIdType>(corners.size()), 1);
	c->GetCellData()->AddArray(decomposition.GetPointer());

	//vtk
	vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	std::stringstream ss;

	//vtu
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	std::stringstream sss;

	//vtm
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	std::stringstream ssss;
	ofstream result;
	//ensight
	vtkSmartPointer<vtkEnSightWriter> wcase = vtkSmartPointer<vtkEnSightWriter>::New();

	switch (output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:
		ss << path << environment->MPIrank << ".vtk";
		writervtk->SetFileName(ss.str().c_str());
		writervtk->SetInputData(c);
		writervtk->Write();
		break;

	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:
		sss << path << environment->MPIrank << ".vtu";
		writervtu->SetFileName(sss.str().c_str());
		writervtu->SetInputData(c);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
		break;

	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:
		ssss << path << environment->MPIrank << ".vtu";
		wvtu->SetFileName(ssss.str().c_str());
		wvtu->SetInputData(c);
		wvtu->SetDataModeToBinary();
		wvtu->Write();

		MPI_Barrier(MPI_COMM_WORLD);
		int size = environment->MPIsize;
		if (environment->MPIrank == 0) {
			result.open("meshC_result.vtm");
			result << "<?xml version=\"1.0\"?>\n";
			result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
			result << " <vtkMultiBlockDataSet>\n";
			for (int i = 0; i < size; i++) {
				result << "  <DataSet index=\"" << i << "\" file=\"meshCorners"
						<< i << ".vtu\">\n  </DataSet>\n";
			}
			result << " </vtkMultiBlockDataSet>\n";
			result << "</VTKFile>\n";
			result.close();
		}
		break;

	case OUTPUT_FORMAT::ENSIGHT_FORMAT:
		vtkSmartPointer<vtkUnstructuredGrid> ugcase = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkMPIController* controller = vtkMPIController::New();
		controller->Initialize();
		//add BlockId
		bool FCD = false;
		if (c->GetCellData()) {
			if (c->GetCellData()->GetArray("BlockId")) {
				FCD = true;
			}
		}
		if (FCD == false) {
			vtkSmartPointer<vtkIntArray> bids = vtkSmartPointer<vtkIntArray>::New();
			bids->SetName("BlockId");
			for (i = 0; i < c->GetNumberOfCells(); i++) {
				bids->InsertNextValue(1);
			}
			c->GetCellData()->SetScalars(bids);
		}
		int blockids[2];
		blockids[0] = 1;
		blockids[1] = 0;
		//write ensight format
		wcase->SetFileName("meshC_result.case");
		wcase->SetNumberOfBlocks(1);
		wcase->SetBlockIDs(blockids);
		wcase->SetTimeStep(0);
		if (environment->MPIrank == 0) {
			int rank = environment->MPIrank;
			controller->Send(c, 0, 1111 + rank);
		}

		if (environment->MPIrank == 0) {
			vtkAppendFilter* app = vtkAppendFilter::New();
			app->AddInputData(c);
			for (int i = 1; i < environment->MPIsize; i++) {
				vtkUnstructuredGrid* h = vtkUnstructuredGrid::New();
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
		break;
	}
	delete[] decomposition_array;
}

void VTK::store(std::vector<std::vector<double> > &displasment)
{
	//storeGeometry(1);
	//storeValues("displacement", 3, displasment, ElementType::NODES);
	return;
	/*
	 //TO DO compresion
	 std::cout << "Compression: " << output->::OUTPUT_COMPRESSION << "\n";
	 //pokus();

	 const std::vector<Element*> &elements = _mesh.elements();
	 const std::vector<eslocal> &_partPtrs = _mesh.getPartition();

	 VTKGrid = vtkUnstructuredGrid::New();

	 size_t n_nodsClust = 0;
	 for (size_t iEl = 0; iEl < elements.size(); iEl++) {
	 n_nodsClust += elements[iEl]->nodes();
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
	 for (size_t j = 0; j < elements[i]->nodes(); j++) {
	 tmp[j] = elements[i]->node(j) + cnt;
	 }
	 VTKGrid->InsertNextCell(elements[i]->vtkCode(), elements[i]->nodes(), &tmp[0]);
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

	 const unsigned int dofs= displasment[0].size() / _mesh.coordinates().localSize(0);

	 double displacement_array[mycounter];

	 counter = 0;



	 for (size_t i = 0; i < displasment.size(); i++) {
	 for (size_t j = 0; j < (displasment[i].size() / dofs); j++) {
	 for(int k=0;k<dofs;k++){
	 displacement_array[dofs * counter + k] = displasment[i][dofs * j + k];
	 }
	 //displacement_array[3 * counter + 1] = displasment[i][3 * j + 1];
	 //displacement_array[3 * counter + 2] = displasment[i][3 * j + 2];

	 counter++;
	 }
	 }

	 std::cout << "NUMBER OF DOFS: " << displasment[0].size() / _mesh.coordinates().localSize(0) << "\n";

	 vtkNew < vtkDoubleArray > displacement;
	 displacement->SetName("displacement");
	 displacement->SetNumberOfComponents(dofs);
	 displacement->SetNumberOfTuples(static_cast<vtkIdType>(counter));
	 VTKGrid->GetPointData()->AddArray(displacement.GetPointer());

	 double* displacementData = displacement_array;
	 vtkIdType numTuples = displacement->GetNumberOfTuples();
	 for (vtkIdType i = 0, counter = 0; i < numTuples; i++, counter++) {
	 //double values[3] = { displacementData[i * 3], displacementData[i * 3 + 1], displacementData[i * 3 + 2] };
	 double values[dofs];
	 for(int k=0;k<dofs;k++){
	 values[k] = {displacementData[i*dofs+k]} ;
	 }

	 displacement->SetTypedTuple(counter, values);
	 }

	 //surface and decimation

	 std::cout << "Decimation: " << output->::OUTPUT_DECIMATION << "\n";
	 vtkSmartPointer < vtkDecimatePro > decimate = vtkSmartPointer
	 < vtkDecimatePro > ::New();
	 vtkAppendFilter* ap = vtkAppendFilter::New();
	 vtkUnstructuredGrid* decimated = vtkUnstructuredGrid::New();
	 vtkGeometryFilter* gf = vtkGeometryFilter::New();
	 vtkDataSetSurfaceFilter* dssf = vtkDataSetSurfaceFilter::New();
	 vtkTriangleFilter* tf = vtkTriangleFilter::New();

	 vtkSmartPointer < vtkContourFilter > cg = vtkSmartPointer < vtkContourFilter
	 > ::New();

	 if (output->::OUTPUT_DECIMATION != 0) {
	 gf->SetInputData(VTKGrid);
	 gf->Update();

	 dssf->SetInputData(gf->GetOutput());
	 dssf->Update();

	 tf->SetInputData(dssf->GetOutput());
	 tf->Update();

	 decimate->SetInputData(tf->GetOutput());
	 decimate->SetTargetReduction(output->::OUTPUT_DECIMATION);
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
	 vtkUnstructuredGrid* ugcase = vtkUnstructuredGrid::New();
	 vtkMPIController* controller = vtkMPIController::New();
	 if(init==false){
	 controller->Initialize();
	 init=true;
	 }
	 int rank = environment->MPIrank;


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

	 switch (output->format) {
	 case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:
	 std::cout << "LEGACY\n";
	 ss << _path << environment->MPIrank<<".vtk";
	 writervtk->SetFileName(ss.str().c_str());
	 if (output->::OUTPUT_DECIMATION == 0) {
	 writervtk->SetInputData(VTKGrid);
	 } else {
	 writervtk->SetInputData(decimated);
	 }
	 if (output->::OUTPUT_COMPRESSION) {
	 // writervtk->SetCompressor(myZlibCompressor);
	 std::cout << "Compression: " << output->::OUTPUT_COMPRESSION
	 << "\n";
	 }

	 writervtk->Write();

	 break;

	 case OUTPUT_FORMAT::VTK_BINARY_FORMAT:
	 std::cout << "BINARY\n";
	 sss << _path << ".vtu";
	 writervtu->SetFileName(sss.str().c_str());
	 if (output->::OUTPUT_DECIMATION == 0) {
	 writervtu->SetInputData(VTKGrid);
	 } else {
	 writervtu->SetInputData(decimated);
	 }
	 writervtu->SetDataModeToBinary();
	 writervtu->Write();
	 break;

	 case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:
	 std::cout << "MULTIBLOCK\n";
	 ssss << _path << ".vtu";
	 wvtu->SetFileName(ssss.str().c_str());
	 if (output->::OUTPUT_DECIMATION == 0) {
	 wvtu->SetInputData(VTKGrid);
	 } else {
	 wvtu->SetInputData(decimated);
	 }
	 wvtu->SetDataModeToBinary();
	 wvtu->Write();

	 MPI_Barrier(MPI_COMM_WORLD);
	 int size = environment->MPIsize;
	 if (environment->MPIrank == 0) {
	 result.open("result.vtm");
	 result << "<?xml version=\"1.0\"?>\n";
	 if (!output->::OUTPUT_COMPRESSION) {
	 result
	 << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	 std::cout << "Compression: "
	 << output->::OUTPUT_COMPRESSION << "\n";
	 } else {
	 result
	 << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
	 std::cout << "Compression: "
	 << output->::OUTPUT_COMPRESSION << "\n";
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

	 case OUTPUT_FORMAT::ENSIGHT_FORMAT:
	 std::cout << "ENSIGHT\n";
	 if (rank != 0) {
	 controller->Send(VTKGrid, 0, 1111 + rank);
	 }

	 //write ensight format
	 wcase->SetFileName("result.case");
	 wcase->SetNumberOfBlocks(1);
	 wcase->SetBlockIDs(blockids);
	 wcase->SetTimeStep(0);
	 if (environment->MPIrank == 0) {
	 vtkAppendFilter* app = vtkAppendFilter::New();
	 app->AddInputData(VTKGrid);
	 for (int i = 1; i < environment->MPIsize; i++) {
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
	 }
	 break;
	 }*/
}

void VTK::gluing(const OutputConfiguration &output, const Mesh &mesh, const Constraints &constraints, const std::string &path, size_t dofs)
{
	 VTK* help=new VTK(output, mesh, path);
	 vtkSmartPointer<vtkPolyData> gp[dofs];
	 vtkSmartPointer<vtkPolyData> dirichlet[dofs];
	 vtkSmartPointer<vtkPoints> po[dofs];
	 vtkSmartPointer<vtkPoints> points[dofs];
	 vtkSmartPointer<vtkPoints> dirpoints[dofs];
	 vtkSmartPointer<vtkCellArray> lines[dofs];
	 vtkSmartPointer<vtkCellArray> ver[dofs];
	 vtkSmartPointer<vtkCellArray> verdir[dofs];
	 for(int i=0;i<dofs;i++){
		 gp[i]=vtkSmartPointer<vtkPolyData>::New();
		 dirichlet[i]=vtkSmartPointer<vtkPolyData>::New();
		 po[i]=vtkSmartPointer<vtkPoints>::New();
		 points[i]=vtkSmartPointer<vtkPoints>::New();
		 lines[i]=vtkSmartPointer<vtkCellArray>::New();
		 ver[i]=vtkSmartPointer<vtkCellArray>::New();
		 verdir[i]=vtkSmartPointer<vtkCellArray>::New();
		 dirpoints[i]=vtkSmartPointer<vtkPoints>::New();
	 }
	 //int V=mesh.coordinates().localSize(0) * constraints.B1.size();
	 std::vector<int> row;
	 std::vector<int> val;

	 MPI_Status status;
	 int matr=0;
	 int l=0;


	 //const Boundaries &sBoundaries = mesh.subdomainBoundaries();
	 //sBoundaries[5] -> vector subdomen
	 //const Boundaries &cBoundaries = mesh.clusterBoundaries();
	 //constraints.B1[0].I_row_indices
	 //constraints.B1[0].J_col_indices
	 //constraints.B1[0].V_values

	 std::vector<int> findex;
	 std::vector<int> fpart;
	 for (int y = 0; y < mesh.parts(); y++) {
		 for(int i=0;i<constraints.B1[y].I_row_indices.size();i++){
			 int h = constraints.B1[y].I_row_indices[i];
			 int index1 = constraints.B1[y].J_col_indices[i]-1;
			 bool dirp = true;
			 for(int p=y;p<mesh.parts();p++){
				 for(int k=0;k<constraints.B1[p].I_row_indices.size();k++){
					 if(constraints.B1[p].I_row_indices[k] == h){
						 if((constraints.B1[y].V_values[i]+constraints.B1[p].V_values[k]) == 0){
							 dirp=false;
							 std::vector<Point> p1(dofs);
							 for(int pd=0;pd<dofs;pd++){
								 p1[pd]=mesh.coordinates().get(index1/dofs,y);
								 p1[pd]=help->shrink(p1[pd],y);
								 vtkIdType pid[1];
								 pid[0] = points[pd]->InsertNextPoint(p1[pd].x,p1[pd].y,p1[pd].z);
								 ver[pd]->InsertNextCell(1,pid);
							 }
							 int index2 = constraints.B1[p].J_col_indices[k]-1;
							 findex.push_back(index2);
							 fpart.push_back(p);
							 Point p2[dofs];
							 for(int pd=0;pd<dofs;pd++){
								 p2[pd] = mesh.coordinates().get(index2/dofs,p);
								 p2[pd] = help->shrink(p2[pd],p);
								 po[pd]->InsertNextPoint(p1[pd].x,p1[pd].y,p1[pd].z);
								 po[pd]->InsertNextPoint(p2[pd].x,p2[pd].y,p2[pd].z);
								 vtkSmartPointer<vtkLine> ls=vtkSmartPointer<vtkLine>::New();
								 ls->GetPointIds()->SetId(0,l);
								 ls->GetPointIds()->SetId(1,l+1);
								 lines[pd]->InsertNextCell(ls);
							 }
							 l += 2;
							 break;
						 }
					 }
				 }
			 }
			 if(dirp == true){
				 for(int i=0;i<findex.size();i++){
					 if(index1 == findex[i] && fpart[i] == y){
						 dirp = false;
						 break;
					 }
				 }
			 }
			 if(dirp == true){
				 row.push_back(h);
				 val.push_back(constraints.B1[y].V_values[i]);
				 matr++;
				 std::vector<Point> p1(dofs);
				 for(int pd=0;pd<dofs;pd++){
					 p1[pd] = mesh.coordinates().get(index1/dofs,y);
					 p1[pd] = help->shrink(p1[pd],y);
					 vtkIdType pid[1];
					 pid[0] = dirpoints[pd]->InsertNextPoint(p1[pd].x,p1[pd].y,p1[pd].z);
					 verdir[pd]->InsertNextCell(1,pid);
				 }
			 }
		 }
	 }



	 //MPI dirichlet control
	 vtkMPIController* controller = vtkMPIController::New();
	 controller->Initialize();
	 vtkSmartPointer<vtkCellArray> vdirMPI[dofs];
	 vtkSmartPointer<vtkPoints> dirpMPI[dofs];
	 for(int i=0;i<dofs;i++){
		vdirMPI[i] = vtkSmartPointer<vtkCellArray>::New();
		dirpMPI[i] = vtkSmartPointer<vtkPoints>::New();
	 }

	 if(environment->MPIsize>1){
		 MPI_Barrier(MPI_COMM_WORLD);

		 vtkSmartPointer<vtkPolyData> PPoints[dofs];
		 for(int i=0;i<dofs;i++){
			 PPoints[i] = vtkSmartPointer<vtkPolyData>::New();
			 PPoints[i]->SetPoints(dirpoints[i]);
		 }

		 int rank = environment->MPIrank;

		 for(int i=0;i<environment->MPIsize;i++){
			 std::vector<int> row2;
			 vtkSmartPointer<vtkPolyData> PPoints2[dofs];
			 for(int pd=0;pd<dofs;pd++){
				 PPoints2[pd]=vtkSmartPointer<vtkPolyData>::New();
			 }
			 std::vector<int> val2;
			 int vs;

			 if(i==rank){
				 row2=row;
				 for(int pd=0;pd<dofs;pd++){
					 PPoints2[pd]->DeepCopy(PPoints[pd]);
				 }
				 val2 = val;
				 vs = row2.size();
			 }
			 MPI_Barrier(MPI_COMM_WORLD);
			 MPI_Bcast(&vs,1,MPI_INT,i,MPI_COMM_WORLD);

			 row2.resize(vs);
			 val2.resize(vs);

			 MPI_Bcast(row2.data(),row2.size(),MPI_INT,i,MPI_COMM_WORLD);
			 MPI_Bcast(val2.data(),val2.size(),MPI_INT,i,MPI_COMM_WORLD);
			 for(int pd=0;pd<dofs;pd++){
				 controller->Broadcast(PPoints2[pd],i);
			 }
			 if(i!=rank){
				 int lg;
				 if(row.size()<row2.size()){
					 lg=row.size();
				 }
				 else{
					 lg=row2.size();
				 }
				 std::vector<int> fi;

				 for(int i=0;i<lg;i++){
					 bool dir = true;
					 for(int y=0;y<lg;y++){
						 if(row[i] == row2[y]){
							 if((val[i] + val2[y]) == 0){
								 dir = false;
								 fi.push_back(row2[y]);
								 for(vtkIdType jed=0;jed<dofs;jed++){
									 double p1[3];
									 double p2[3];
									 PPoints[jed]->GetPoint(i,p1);
									 PPoints2[jed]->GetPoint(y,p2);
									 po[jed]->InsertNextPoint(p1[0],p1[1],p1[2]);
									 po[jed]->InsertNextPoint(p2[0],p2[1],p2[2]);

									 vtkSmartPointer<vtkLine> ls=vtkSmartPointer<vtkLine>::New();
									 ls->GetPointIds()->SetId(0,l);
									 ls->GetPointIds()->SetId(1,l+1);
									 lines[jed]->InsertNextCell(ls);
								 }
								 l += 2;
								 break;
							 }
						 }
					 }
					 if(dir == true){
						 for(int n=0;n<fi.size();n++){
							 if(row[i] == fi[n]){
								 dir = false;
								 break;
							 }
						 }
						 if(dir == true){
							 for(size_t t=0;t<constraints.B1clustersMap.size();t++){
								 if((constraints.B1clustersMap[t][0]+1)==row[i] && constraints.B1clustersMap[t].size()==2){
									 for(vtkIdType jed=0;jed<dofs;jed++){
										 double p1[3];
										 PPoints[jed]->GetPoint(i,p1);
										 vtkIdType pid[1];
										 pid[0] = dirpMPI[jed]->InsertNextPoint(p1[0],p1[1],p1[2]);
										 vdirMPI[jed]->InsertNextCell(1,pid);
									 }
									 break;
								 }
							 }
						 }
					 }
				 }
			 }
			 MPI_Barrier(MPI_COMM_WORLD);
		 }
	 }
	 vtkSmartPointer<vtkAppendFilter> ap=vtkSmartPointer<vtkAppendFilter>::New();
	 vtkSmartPointer<vtkAppendFilter> app=vtkSmartPointer<vtkAppendFilter>::New();

	 for(size_t i=0;i<dofs;i++){
		 if(environment->MPIsize==1){
			 dirichlet[i]->SetPoints(dirpoints[i]);
			 dirichlet[i]->SetVerts(verdir[i]);
		 }else{
			 dirichlet[i]->SetPoints(dirpMPI[i]);
			 dirichlet[i]->SetVerts(vdirMPI[i]);
		 }

		 app->AddInputData(dirichlet[i]);

		 gp[i]->SetPoints(po[i]);
		 gp[i]->SetLines(lines[i]);

		 ap->AddInputData(gp[i]);

	 }
	 app->Update();
	 ap->Update();

	 stringstream ss;
	 switch (output.format) {
	 	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
	 		vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	 		ss<<path<<environment->MPIrank<<".vtk";
	 		writervtk->SetFileName(ss.str().c_str());
	 		writervtk->SetInputData(ap->GetOutput());
	 		writervtk->Write();
	 		vtkSmartPointer<vtkGenericDataObjectWriter> dw=vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	 		stringstream sd;
	 		sd<<"dirichlet"<<environment->MPIrank<<".vtk";
	 		dw->SetFileName(sd.str().c_str());
	 		dw->SetInputData(app->GetOutput());
	 		dw->Write();
	 	}
	 	break;
	 	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
	 		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	 		ss<<path<<environment->MPIrank<<".vtu";
	 		writervtu->SetFileName(ss.str().c_str());
	 		writervtu->SetInputData(ap->GetOutput());
	 		writervtu->SetDataModeToBinary();
	 		writervtu->Write();
	 		vtkSmartPointer<vtkXMLUnstructuredGridWriter> dw1=vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	 		stringstream sd;
	 		sd<<"dirichlet"<<environment->MPIrank<<".vtu";
	 		dw1->SetFileName(sd.str().c_str());
	 		dw1->SetInputData(app->GetOutput());
	 		dw1->Write();
	 	}
	 	break;
	 	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
	 		vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	 		std::ofstream result;
	 		ss<<path<<environment->MPIrank<<".vtu";
	 		wvtu->SetFileName(ss.str().c_str());
	 		wvtu->SetInputData(ap->GetOutput());
	 		wvtu->SetDataModeToBinary();
	 		wvtu->Write();
	 		vtkSmartPointer<vtkXMLUnstructuredGridWriter> dw2=vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	 		stringstream sd;
	 		sd<<"dirichlet"<<environment->MPIrank<<".vtu";
	 		dw2->SetFileName(sd.str().c_str());
	 		dw2->SetInputData(app->GetOutput());
	 		dw2->Write();

	 		int size = environment->MPIsize;
	 		if (environment->MPIrank == 0) {
	 			stringstream sss;
	 			sss << path << ".vtm";
	 			result.open(sss.str().c_str());
	 			result << "<?xml version=\"1.0\"?>\n";
	 			if (!output.compression) {
	 				result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	 			} else {
	 				result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
	 			}
	 			result << " <vtkMultiBlockDataSet>\n";
	 			int vtmi=0;
	 			for (int i = 0; i < size; i++) {
	 				result << "  <DataSet index=\"" << vtmi << "\" file=\"" << path << i << ".vtu\">\n  </DataSet>\n";
	 				result << "  <DataSet index=\"" << vtmi+1 << "\" file=\"" << "dirichlet" << i << ".vtu\">\n  </DataSet>\n";
	 				vtmi+=2;
	 			}
	 			result << " </vtkMultiBlockDataSet>\n";
	 			result << "</VTKFile>\n";
	 			result.close();
	 		}
	 	}
	 	break;
	 	case OUTPUT_FORMAT::ENSIGHT_FORMAT: {
	 		vtkSmartPointer<vtkEnSightWriter> wcase = vtkSmartPointer<vtkEnSightWriter>::New();
	 		vtkSmartPointer<vtkUnstructuredGrid> ugcase = vtkSmartPointer<vtkUnstructuredGrid>::New();
	 		vtkSmartPointer<vtkUnstructuredGrid> VTKGrid1 = vtkSmartPointer<vtkUnstructuredGrid>::New();
	 		vtkSmartPointer<vtkUnstructuredGrid> VTKGrid2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
	 		VTKGrid1->ShallowCopy(ap->GetOutput());
	 		VTKGrid2->ShallowCopy(app->GetOutput());
	 		bool FCD = false;
	 		if (VTKGrid1->GetCellData()) {
	 			if (VTKGrid1->GetCellData()->GetArray("BlockId")) {
	 				FCD = true;
	 			}
	 		}
	 		if (FCD == false) {
	 			vtkSmartPointer<vtkIntArray> bids = vtkSmartPointer<vtkIntArray>::New();
	 			bids->SetName("BlockId");
	 			for (int i = 0; i < VTKGrid1->GetNumberOfCells(); i++) {
	 				bids->InsertNextValue(1);
	 			}
	 			VTKGrid1->GetCellData()->SetScalars(bids);
	 		}
	 		FCD = false;
	 		if (VTKGrid2->GetCellData()) {
	 			if (VTKGrid2->GetCellData()->GetArray("BlockId")) {
	 				FCD = true;
	 			}
	 		}
	 		if (FCD == false) {
	 			vtkSmartPointer<vtkIntArray> bids = vtkSmartPointer<vtkIntArray>::New();
	 			bids->SetName("BlockId");
	 			for (int i = 0; i < VTKGrid2->GetNumberOfCells(); i++) {
	 				bids->InsertNextValue(1);
	 			}
	 			VTKGrid2->GetCellData()->SetScalars(bids);
	 		}
	 		int blockids[2];
	 		blockids[0] = 1;
	 		blockids[1] = 0;

	 		if (environment->MPIrank != 0) {
	 			controller->Send(VTKGrid1, 0, 1111 + environment->MPIrank);
	 			controller->Send(VTKGrid2, 0, 2222 + environment->MPIrank);
	 		}

	 		if (environment->MPIrank == 0) {
	 			stringstream sc;
	 			sc<<path<<".case";
	 			wcase->SetFileName(sc.str().c_str());
	 			wcase->SetNumberOfBlocks(1);
	 			wcase->SetBlockIDs(blockids);
	 			wcase->SetTimeStep(0);
	 			vtkSmartPointer<vtkAppendFilter> appp = vtkSmartPointer<vtkAppendFilter>::New();
	 			appp->AddInputData(VTKGrid1);
	 			appp->AddInputData(VTKGrid2);
	 			for (int i = 1; i < environment->MPIsize; i++) {
	 				vtkSmartPointer<vtkUnstructuredGrid> h = vtkSmartPointer<vtkUnstructuredGrid>::New();
	 				vtkSmartPointer<vtkUnstructuredGrid> h2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
	 				controller->Receive(h, i, 1111 + i);
	 				controller->Receive(h2, i, 2222 + i);
	 				appp->AddInputData(h);
	 				appp->AddInputData(h2);
	 			}
	 			appp->Update();
	 			//std::cout<<app->GetOutput()->GetNumberOfPoints()<<std::endl;
	 			vtkSmartPointer<vtkGeometryFilter> gf = vtkSmartPointer<vtkGeometryFilter>::New();
	 			gf->SetInputData(appp->GetOutput());
	 			gf->Update();
	 			vtkSmartPointer<vtkCleanPolyData> cpd = vtkSmartPointer<vtkCleanPolyData>::New();
	 			cpd->SetInputData(gf->GetOutput());
	 			cpd->Update();
	 			vtkSmartPointer<vtkAppendFilter> apc = vtkSmartPointer<vtkAppendFilter>::New();
	 			apc->SetInputData(cpd->GetOutput());
	 			apc->Update();
	 			//std::cout<<apc->GetOutput()->GetNumberOfPoints()<<std::endl;
	 			ugcase->ShallowCopy(apc->GetOutput());
	 			wcase->SetInputData(ugcase);
	 			wcase->Write();
	 			wcase->WriteCaseFile(1);
	 		}
	 	}
	 	break;
	 }

	 delete help;
	 controller->Delete();
}

void VTK::finalize(){
	vtkSmartPointer<vtkUnstructuredGrid> VTKGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	stringstream ss;
	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
		vtkSmartPointer<vtkGenericDataObjectReader> rvtk = vtkSmartPointer<vtkGenericDataObjectReader>::New();
		ss<<_path<<environment->MPIrank<<".vtk";
		rvtk->SetFileName(ss.str().c_str());
		rvtk->Update();
		VTKGrid->ShallowCopy(rvtk->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rvtu = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rvtu->SetFileName(ss.str().c_str());
		rvtu->Update();
		VTKGrid->ShallowCopy(rvtu->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rvtm = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rvtm->SetFileName(ss.str().c_str());
		rvtm->Update();
		VTKGrid->ShallowCopy(rvtm->GetOutput());
	}
	break;
	case OUTPUT_FORMAT::ENSIGHT_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridReader> rcase = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
		ss<<_path<<environment->MPIrank<<".vtu";
		rcase->SetFileName(ss.str().c_str());
		rcase->Update();
		VTKGrid->ShallowCopy(rcase->GetOutput());
	}
	break;
	}

	vtkSmartPointer < vtkDecimatePro > decimate = vtkSmartPointer < vtkDecimatePro > ::New();
	vtkSmartPointer<vtkAppendFilter> ap = vtkSmartPointer<vtkAppendFilter>::New();
	vtkSmartPointer<vtkGeometryFilter> gf = vtkSmartPointer<vtkGeometryFilter>::New();
	vtkSmartPointer<vtkDataSetSurfaceFilter> dssf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	vtkSmartPointer<vtkTriangleFilter> tf = vtkSmartPointer<vtkTriangleFilter>::New();

	if (_output.decimation != 0) {
		 gf->SetInputData(VTKGrid);
		 gf->Update();

		 dssf->SetInputData(gf->GetOutput());
		 dssf->Update();

		 tf->SetInputData(dssf->GetOutput());
		 tf->Update();

		 decimate->SetInputData(tf->GetOutput());
		 decimate->SetTargetReduction(_output.decimation);
		 decimate->Update();

		 ap->AddInputData(decimate->GetOutput());
		 ap->Update();

		 VTKGrid->ShallowCopy(ap->GetOutput());
	}

	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT:{
		vtkSmartPointer<vtkGenericDataObjectWriter> writervtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
		writervtk->SetFileName(ss.str().c_str());
		writervtk->SetInputData(VTKGrid);
		writervtk->Write();
	}
	break;
	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writervtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writervtu->SetFileName(ss.str().c_str());
		writervtu->SetInputData(VTKGrid);
		writervtu->SetDataModeToBinary();
		writervtu->Write();
	}
	break;
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT:{
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> wvtu = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		wvtu->SetFileName(ss.str().c_str());
		wvtu->SetInputData(VTKGrid);
		wvtu->SetDataModeToBinary();
		wvtu->Write();
	}
	break;
	case OUTPUT_FORMAT::ENSIGHT_FORMAT: {
		vtkMPIController* controller=vtkMPIController::New();
		controller->Initialize();
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

		if (environment->MPIrank != 0) {
			controller->Send(VTKGrid, 0, 1111 + environment->MPIrank);
		}


		if (environment->MPIrank == 0) {
			wcase->SetFileName(ss.str().c_str());
			wcase->SetNumberOfBlocks(1);
			wcase->SetBlockIDs(blockids);
			wcase->SetTimeStep(0);
			vtkSmartPointer<vtkAppendFilter> app = vtkSmartPointer<vtkAppendFilter>::New();
			app->AddInputData(VTKGrid);
			for (int i = 1; i < environment->MPIsize; i++) {
				vtkSmartPointer<vtkUnstructuredGrid> h = vtkSmartPointer<vtkUnstructuredGrid>::New();
				controller->Receive(h, i, 1111 + i);
				app->AddInputData(h);
			}
			app->Update();
			//std::cout<<app->GetOutput()->GetNumberOfPoints()<<std::endl;
			vtkSmartPointer<vtkGeometryFilter> gf = vtkSmartPointer<vtkGeometryFilter>::New();
			gf->SetInputData(app->GetOutput());
			gf->Update();
			vtkSmartPointer<vtkCleanPolyData> cpd = vtkSmartPointer<vtkCleanPolyData>::New();
			cpd->SetInputData(gf->GetOutput());
			cpd->Update();
			vtkSmartPointer<vtkAppendFilter> apc = vtkSmartPointer<vtkAppendFilter>::New();
			apc->SetInputData(cpd->GetOutput());
			apc->Update();
			//std::cout<<apc->GetOutput()->GetNumberOfPoints()<<std::endl;
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
