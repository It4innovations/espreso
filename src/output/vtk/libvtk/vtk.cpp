
#include "vtkNew.h"
#include "vtkSmartPointer.h"

#include "vtkIntArray.h"
#include "vtkDoubleArray.h"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkUnstructuredGrid.h"

#include "vtkGeometryFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkTriangleFilter.h"
#include "vtkDecimatePro.h"
#include "vtkAppendFilter.h"

#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkGenericDataObjectWriter.h"
#include "vtkEnSightWriter.h"

#include "../vtk.h"

using namespace espreso::store;

VTK::VTK(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: Store(output, mesh, path), _lastData(ElementType::ELEMENTS)
{
	computeCenters();
	VTKGrid = vtkUnstructuredGrid::New();
}

VTK::VTK(const OutputConfiguration &output, const Mesh &mesh)
: Store(output, mesh, ""), _lastData(ElementType::ELEMENTS)
{
	computeCenters();
	VTKGrid = vtkUnstructuredGrid::New();
}

VTK::~VTK()
{
	for (size_t i = 0; i < VTKDataArrays.size(); i++) {
		delete[] VTKDataArrays[i];
	}
	VTKGrid->Delete();
}

void VTK::coordinates()
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t d = 0; d < _mesh.parts(); d++) {
		for (size_t i = 0; i < _mesh.coordinates().localSize(d); i++) {
			espreso::Point p = shrink(_mesh.coordinates().get(i, d), d);
			points->InsertNextPoint(p.x, p.y, p.z);
		}
	}
	VTKGrid->SetPoints(points);
}

static void storeDecomposition(vtkUnstructuredGrid *VTKGrid, std::vector<void*> &VTKDataArrays, eslocal *decomposition, size_t size)
{
	vtkNew<vtkIntArray> vtkDecomposition;
	vtkDecomposition->SetName("decomposition");
	vtkDecomposition->SetNumberOfComponents(1);
	vtkDecomposition->SetArray(decomposition, static_cast<vtkIdType>(size), 1);
	VTKGrid->GetCellData()->AddArray(vtkDecomposition.GetPointer());
	VTKDataArrays.push_back(decomposition);
}

void VTK::nodes(const std::vector<Element*> &nodes)
{
	size_t nSize = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		nSize += nodes[i]->domains().size();
	}

	eslocal *decomposition = new eslocal[nSize];
	VTKGrid->Allocate(static_cast<vtkIdType>(nSize));

	vtkIdType node = 0;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t i = 0, c = 0; i < nodes.size(); i++) {
		for (size_t d = 0; d < nodes[i]->domains().size(); d++, node++) {
			espreso::Point p = shrink(_mesh.coordinates()[nodes[i]->node(0)], nodes[i]->domains()[d]);
			points->InsertNextPoint(p.x, p.y, p.z);
			decomposition[c++] = nodes[i]->domains()[d];
			VTKGrid->InsertNextCell(nodes[i]->vtkCode(), 1, &node);
		}
	}
	VTKGrid->SetPoints(points);

	storeDecomposition(VTKGrid, VTKDataArrays, decomposition, nSize);
}

void VTK::nodes(const std::vector<std::vector<eslocal> > &nodes)
{
	size_t nSize = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		nSize += nodes[i].size();
	}

	eslocal *decomposition = new eslocal[nSize];
	VTKGrid->Allocate(static_cast<vtkIdType>(nSize));

	vtkIdType node = 0;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t part = 0, c = 0; part < nodes.size(); part++) {
		for (size_t n = 0; n < nodes[part].size(); n++, node++) {
			espreso::Point p = shrink(_mesh.coordinates()[nodes[part][n]], part);
			points->InsertNextPoint(p.x, p.y, p.z);
			decomposition[c++] = part;
			VTKGrid->InsertNextCell(NodeVTKCode, 1, &node);
		}
	}
	VTKGrid->SetPoints(points);

	storeDecomposition(VTKGrid, VTKDataArrays, decomposition, nSize);
}

void VTK::cells(ElementType eType)
{
	std::vector<espreso::Element*> elements;

	switch (eType) {
	case espreso::store::Store::ElementType::ELEMENTS:
		elements.insert(elements.end(), _mesh.elements().begin(), _mesh.elements().end());
		break;
	case espreso::store::Store::ElementType::FACES:
		elements.insert(elements.end(), _mesh.faces().begin(), _mesh.faces().end());
		std::sort(elements.begin(), elements.end(), [] (const espreso::Element* e1, const espreso::Element *e2) { return e1->domains() < e2->domains(); });
		break;
	case espreso::store::Store::ElementType::EDGES:
		elements.insert(elements.end(), _mesh.edges().begin(), _mesh.edges().end());
		std::sort(elements.begin(), elements.end(), [] (const espreso::Element* e1, const espreso::Element *e2) { return e1->domains() < e2->domains(); });
		break;
	default:
		break;
	}

	size_t nSize = 0, eSize = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		nSize += elements[i]->domains().size() * elements[i]->nodes();
		eSize += elements[i]->domains().size();
	}

	std::vector<size_t> offset = { 0 };
	for (size_t p = 1; p < _mesh.parts(); p++) {
		offset.push_back(offset[p - 1] + _mesh.coordinates().localSize(p - 1));
	}

	std::vector<vtkIdType> nodes(20);
	eslocal *decomposition = new eslocal[eSize];

	VTKGrid->Allocate(static_cast<vtkIdType>(nSize));

	for (size_t i = 0, c = 0; i < elements.size(); i++) {
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			nodes.clear();
			for (size_t n = 0; n < elements[i]->nodes(); n++) {
				nodes.push_back(_mesh.coordinates().localIndex(elements[i]->node(n), elements[i]->domains()[d]) + offset[elements[i]->domains()[d]]);
			}
			decomposition[c++] = elements[i]->domains()[d];
			VTKGrid->InsertNextCell(elements[i]->vtkCode(), elements[i]->nodes(), nodes.data());
		}
	}

	storeDecomposition(VTKGrid, VTKDataArrays, decomposition, eSize);
}

void VTK::lambdas(const std::vector<std::vector<eslocal> > &nodes, std::function<Point(const Point&, size_t, size_t, bool)> shrink)
{
	size_t nSize = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		nSize += nodes[i].size();
	}

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t part = 0; part < nodes.size(); part++) {
		for (size_t n = 0; n < nodes[part].size(); n++) {
			espreso::Point p = shrink(_mesh.coordinates()[nodes[part][n]], part, n, true);
			points->InsertNextPoint(p.x, p.y, p.z);
			p = shrink(_mesh.coordinates()[nodes[part][n]], part, n, false);
			points->InsertNextPoint(p.x, p.y, p.z);
		}
	}
	VTKGrid->SetPoints(points);

	std::vector<vtkIdType> line = { 0, 1 };
	eslocal *decomposition = new eslocal[nSize];

	VTKGrid->Allocate(static_cast<vtkIdType>(2 * nSize));

	for (size_t part = 0, c = 0; part < nodes.size(); part++) {
		for (size_t n = 0; n < nodes[part].size(); n++) {
			decomposition[c++] = part;
			VTKGrid->InsertNextCell(Line2VTKCode, 2, line.data());
			line[0] += 2;
			line[1] += 2;
		}
	}

	storeDecomposition(VTKGrid, VTKDataArrays, decomposition, nSize);
}

template <typename TVTKArray, typename TType>
static void storeData(vtkUnstructuredGrid *VTKGrid, std::vector<void*> &VTKDataArrays, const std::string &name, size_t dimension, const std::vector<std::vector<TType> > &values, espreso::store::Store::ElementType eType)
{
	size_t size = 0;
	for (size_t i = 0; i < values.size(); i++) {
		size += values[i].size();
	}

	TType *data = new TType[size];
	for (size_t i = 0, offset = 0; i < values.size(); offset += values[i++].size()) {
		memcpy(data + offset, values[i].data(), values[i].size() * sizeof(TType));
	}

	vtkNew<TVTKArray> vtkArray;
	vtkArray->SetName(name.c_str());
	vtkArray->SetNumberOfComponents(dimension);
	vtkArray->SetArray(data, static_cast<vtkIdType>(size), 1);

	switch (eType) {
	case espreso::store::Store::ElementType::NODES:
		VTKGrid->GetPointData()->AddArray(vtkArray.GetPointer());
		if (VTKGrid->GetPointData()->GetNumberOfArrays() == 1) {
			VTKGrid->GetPointData()->SetActiveScalars(name.c_str());
		}
		break;
	case espreso::store::Store::ElementType::ELEMENTS:
		VTKGrid->GetCellData()->AddArray(vtkArray.GetPointer());
		if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0) {
			VTKGrid->GetCellData()->SetActiveScalars(name.c_str());
		}
		break;
	}
	VTKDataArrays.push_back(data);
}

void VTK::data(const std::string &name, size_t dimension, const std::vector<std::vector<eslocal> > &values, espreso::store::Store::ElementType eType)
{
	storeData<vtkIntArray, eslocal>(VTKGrid, VTKDataArrays, name, dimension, values, eType);
}

void VTK::data(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, espreso::store::Store::ElementType eType)
{
	storeData<vtkDoubleArray, double>(VTKGrid, VTKDataArrays, name, dimension, values, eType);
}

void VTK::finalize()
{
	std::stringstream name;
	name << _path << espreso::environment->MPIrank;

	if (_output.decimation) {
		vtkSmartPointer<vtkGeometryFilter>       geometry  = vtkSmartPointer<vtkGeometryFilter>::New();
		vtkSmartPointer<vtkDataSetSurfaceFilter> surface   = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
		vtkSmartPointer<vtkTriangleFilter>       triangles = vtkSmartPointer<vtkTriangleFilter>::New();
		vtkSmartPointer<vtkDecimatePro>          decimated = vtkSmartPointer<vtkDecimatePro>::New();
		vtkSmartPointer<vtkAppendFilter>         result    = vtkSmartPointer<vtkAppendFilter>::New();

		geometry->SetInputData(VTKGrid);
		geometry->Update();

		surface->SetInputConnection(geometry->GetOutputPort());
		surface->Update();

		triangles->SetInputConnection(surface->GetOutputPort());
		triangles->Update();

		decimated->SetInputConnection(triangles->GetOutputPort());
		decimated->SetTargetReduction(_output.decimation);
		decimated->Update();

		result->AddInputConnection(decimated->GetOutputPort());
		result->Update();

		VTKGrid->ShallowCopy(result->GetOutput());
	}

	switch (_output.format) {
	case OUTPUT_FORMAT::VTK_LEGACY_FORMAT: {
		vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
		writer->SetFileName((name.str() + ".vtk").c_str());
		writer->SetInputData(VTKGrid);
		writer->Write();
	} break;

	case OUTPUT_FORMAT::VTK_BINARY_FORMAT:
	case OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT: {
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName((name.str() + ".vtu").c_str());
		writer->SetInputData(VTKGrid);
		writer->SetDataModeToBinary();
		if (_output.compression) {
			writer->SetCompressorTypeToZLib();
		} else {
			writer->SetCompressorTypeToNone();
		}
		writer->Write();

		if (_output.format == OUTPUT_FORMAT::VTK_BINARY_FORMAT || espreso::environment->MPIrank) {
			break;
		}

		std::ofstream result(_path + ".vtm");

		result << "<?xml version=\"1.0\"?>\n";
		if (_output.compression) {
			result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
		} else {
			result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
		}
		result << " <vtkMultiBlockDataSet>\n";
		for (int i = 0; i < espreso::environment->MPIsize; i++) {
			result << "  <DataSet index=\"" << i << "\" file=\"" << _path << i << ".vtu\"> </DataSet>\n";
		}
		result << " </vtkMultiBlockDataSet>\n";
		result << "</VTKFile>\n";
		result.close();
	} break;

	case OUTPUT_FORMAT::ENSIGHT_FORMAT: {
		if (!VTKGrid->GetCellData()->GetArray("BlockId")) {
			// EnSight needs Block ID for each cell
			std::vector<std::vector<eslocal> > blockID(_mesh.parts());
			for (size_t p = 0; p < _mesh.parts(); p++) {
				blockID[p].resize(_mesh.getPartition()[p + 1] - _mesh.getPartition()[p], 1);
			}
			storeData<vtkIntArray, eslocal>(VTKGrid, VTKDataArrays, "BlockId", 1, blockID, ElementType::ELEMENTS);
		}

		vtkSmartPointer<vtkEnSightWriter> writer = vtkSmartPointer<vtkEnSightWriter>::New();
		writer->SetFileName(name.str().c_str());
		writer->SetInputData(VTKGrid);
		writer->Write();
		writer->WriteCaseFile(1);
	} break;
	}
}
