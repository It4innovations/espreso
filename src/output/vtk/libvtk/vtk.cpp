
#include "vtkSmartPointer.h"
#include "vtkNew.h"

#include "vtkIntArray.h"
#include "vtkDoubleArray.h"

#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

#include "vtkGenericDataObjectWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include "vtkGenericDataObjectReader.h"
#include "vtkXMLUnstructuredGridReader.h"

#include "../vtk.h"

using namespace espreso::store;

VTK::VTK(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: Store(output, mesh, path), _lastData(ElementType::ELEMENTS)
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

static void points(vtkUnstructuredGrid *VTKGrid, const espreso::Mesh &mesh, std::function<espreso::Point(const espreso::Point&, size_t)> shrink)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t d = 0; d < mesh.parts(); d++) {
		for (size_t i = 0; i < mesh.coordinates().localSize(d); i++) {
			espreso::Point p = shrink(mesh.coordinates().get(i, d), d);
			points->InsertNextPoint(p.x, p.y, p.z);
		}
	}
	VTKGrid->SetPoints(points);
}

static void points(vtkUnstructuredGrid *VTKGrid, std::vector<void*> &VTKDataArrays, const espreso::Mesh &mesh, const std::vector<espreso::Element*> &nodes, std::function<espreso::Point(const espreso::Point&, size_t)> shrink)
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
			espreso::Point p = shrink(mesh.coordinates()[nodes[i]->node(0)], nodes[i]->domains()[d]);
			points->InsertNextPoint(p.x, p.y, p.z);
			decomposition[c++] = nodes[i]->domains()[d];
			VTKGrid->InsertNextCell(nodes[i]->vtkCode(), 1, &node);
		}
	}
	VTKGrid->SetPoints(points);

	vtkNew<vtkIntArray> vtkDecomposition;
	vtkDecomposition->SetName("decomposition");
	vtkDecomposition->SetNumberOfComponents(1);
	vtkDecomposition->SetArray(decomposition, static_cast<vtkIdType>(nSize), 1);
	VTKGrid->GetCellData()->AddArray(vtkDecomposition.GetPointer());
	VTKDataArrays.push_back(decomposition);
}

static void points(vtkUnstructuredGrid *VTKGrid, std::vector<void*> &VTKDataArrays, const espreso::Mesh &mesh, const std::vector<std::vector<eslocal> > &nodes, std::function<espreso::Point(const espreso::Point&, size_t)> shrink)
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
			espreso::Point p = shrink(mesh.coordinates()[nodes[part][n]], part);
			points->InsertNextPoint(p.x, p.y, p.z);
			decomposition[c++] = part;
			VTKGrid->InsertNextCell(NodeVTKCode, 1, &node);
		}
	}
	VTKGrid->SetPoints(points);

	vtkNew<vtkIntArray> vtkDecomposition;
	vtkDecomposition->SetName("decomposition");
	vtkDecomposition->SetNumberOfComponents(1);
	vtkDecomposition->SetArray(decomposition, static_cast<vtkIdType>(nSize), 1);
	VTKGrid->GetCellData()->AddArray(vtkDecomposition.GetPointer());
	VTKDataArrays.push_back(decomposition);
}

static void cells(vtkUnstructuredGrid *VTKGrid, std::vector<void*> &VTKDataArrays, const espreso::Mesh &mesh, espreso::store::Store::ElementType eType)
{
	std::vector<espreso::Element*> elements;

	switch (eType) {
	case espreso::store::Store::ElementType::ELEMENTS:
		elements.insert(elements.end(), mesh.elements().begin(), mesh.elements().end());
		break;
	case espreso::store::Store::ElementType::FACES:
		elements.insert(elements.end(), mesh.faces().begin(), mesh.faces().end());
		std::sort(elements.begin(), elements.end(), [] (const espreso::Element* e1, const espreso::Element *e2) { return e1->domains() < e2->domains(); });
		break;
	case espreso::store::Store::ElementType::EDGES:
		elements.insert(elements.end(), mesh.edges().begin(), mesh.edges().end());
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
	for (size_t p = 1; p < mesh.parts(); p++) {
		offset.push_back(offset[p - 1] + mesh.coordinates().localSize(p - 1));
	}

	std::vector<vtkIdType> nodes(20);
	eslocal *decomposition = new eslocal[eSize];

	VTKGrid->Allocate(static_cast<vtkIdType>(nSize));

	for (size_t i = 0, c = 0; i < elements.size(); i++) {
		for (size_t d = 0; d < elements[i]->domains().size(); d++) {
			nodes.clear();
			for (size_t n = 0; n < elements[i]->nodes(); n++) {
				nodes.push_back(mesh.coordinates().localIndex(elements[i]->node(n), elements[i]->domains()[d]) + offset[elements[i]->domains()[d]]);
			}
			decomposition[c++] = elements[i]->domains()[d];
			VTKGrid->InsertNextCell(elements[i]->vtkCode(), elements[i]->nodes(), nodes.data());
		}
	}

	vtkNew<vtkIntArray> vtkDecomposition;
	vtkDecomposition->SetName("decomposition");
	vtkDecomposition->SetNumberOfComponents(1);
	vtkDecomposition->SetArray(decomposition, static_cast<vtkIdType>(eSize), 1);
	VTKGrid->GetCellData()->AddArray(vtkDecomposition.GetPointer());
	VTKDataArrays.push_back(decomposition);
}

static void lambdas(vtkUnstructuredGrid *VTKGrid, std::vector<void*> &VTKDataArrays, const espreso::Mesh &mesh, const std::vector<std::vector<eslocal> > &nodes, std::function<espreso::Point(const espreso::Point&, size_t, size_t, bool)> shrink)
{
	size_t nSize = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		nSize += nodes[i].size();
	}

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (size_t part = 0; part < nodes.size(); part++) {
		for (size_t n = 0; n < nodes[part].size(); n++) {
			espreso::Point p = shrink(mesh.coordinates()[nodes[part][n]], part, n, true);
			points->InsertNextPoint(p.x, p.y, p.z);
			p = shrink(mesh.coordinates()[nodes[part][n]], part, n, false);
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

	vtkNew<vtkIntArray> vtkDecomposition;
	vtkDecomposition->SetName("decomposition");
	vtkDecomposition->SetNumberOfComponents(1);
	vtkDecomposition->SetArray(decomposition, static_cast<vtkIdType>(nSize), 1);
	VTKGrid->GetCellData()->AddArray(vtkDecomposition.GetPointer());
	VTKDataArrays.push_back(decomposition);
}

template <typename TType, typename TVTKType>
static void data(vtkUnstructuredGrid *VTKGrid, std::vector<void*> &VTKDataArrays, const std::string &name, size_t dimension, const std::vector<std::vector<TType> > &values, espreso::store::Store::ElementType eType)
{
	size_t size = 0;
	for (size_t i = 0; i < values.size(); i++) {
		size += values[i].size();
	}

	TType *data = new TType[size];
	for (size_t i = 0, offset = 0; i < values.size(); offset += values[i++].size()) {
		memcpy(data + offset, values[i].data(), values[i].size() * sizeof(TType));
	}

	vtkNew<TVTKType> vtkArray;
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
		break;
	}
	VTKDataArrays.push_back(data);
}

static void save(vtkUnstructuredGrid *VTKGrid, const std::string &path, const espreso::OutputConfiguration &configuration)
{
	std::stringstream name;
	name << path << espreso::environment->MPIrank;

	switch (configuration.format) {
	case espreso::OUTPUT_FORMAT::VTK_LEGACY_FORMAT: {
		vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
		writer->SetFileName((name.str() + ".vtk").c_str());
		writer->SetInputData(VTKGrid);
		writer->Write();
	} break;

	case espreso::OUTPUT_FORMAT::VTK_BINARY_FORMAT: {
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName((name.str() + ".vtu").c_str());
		writer->SetInputData(VTKGrid);
		writer->SetDataModeToBinary();
		writer->Write();
	} break;

	case espreso::OUTPUT_FORMAT::VTK_MULTIBLOCK_FORMAT: {
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

		writer->SetFileName((name.str() + ".vtu").c_str());
		writer->SetInputData(VTKGrid);
		writer->SetDataModeToBinary();
		writer->Write();

		if (espreso::environment->MPIrank) {
			break;
		}
		std::ofstream result(path + ".vtm");

		result << "<?xml version=\"1.0\"?>\n";
		if (!configuration.compression) {
			result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
		} else {
			result << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n";
		}
		result << " <vtkMultiBlockDataSet>\n";
		for (int i = 0; i < espreso::environment->MPIsize; i++) {
			result << "  <DataSet index=\"" << i << "\" file=\"" << path << i << ".vtu\"> </DataSet>\n";
		}
		result << " </vtkMultiBlockDataSet>\n";
		result << "</VTKFile>\n";
		result.close();
	} break;

	case espreso::OUTPUT_FORMAT::ENSIGHT_FORMAT: {
		ESINFO(espreso::GLOBAL_ERROR) << "Implement ENSIGHT";
	} break;
	}
}

void VTK::storeGeometry(size_t timeStep)
{
	points(VTKGrid, _mesh, [&] (const Point &point, size_t part) { return shrink(point, part); });
	cells(VTKGrid, VTKDataArrays, _mesh, ElementType::ELEMENTS);
}

void VTK::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{
	std::vector<std::vector<int> > selection(_mesh.parts());
	std::vector<std::vector<double> > values(_mesh.parts());

	switch (eType) {
	case ElementType::ELEMENTS:
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
		data<double, vtkDoubleArray>(VTKGrid, VTKDataArrays, name + "FixedValue", properties.size(), values, eType);
		break;
	case ElementType::FACES:
		ESINFO(GLOBAL_ERROR) << "Implement store properties for faces";
		break;
	case ElementType::EDGES:
		ESINFO(GLOBAL_ERROR) << "Implement store properties for edges";
		break;
	case ElementType::NODES:
		for (size_t n = 0; n < _mesh.nodes().size(); n++) {
			const Element *node = _mesh.nodes()[n];
			for (size_t p = 0; p < properties.size(); p++) {
				for (size_t d = 0; d < node->domains().size(); d++) {
					size_t domain = node->domains()[d];
					selection[domain].push_back(node->settings().isSet(properties[p]) ? 1 : 0);
					values[domain].push_back(node->settings(properties[p]).back()->evaluate(n));
				}
			}
		}
		data<int, vtkIntArray>(VTKGrid, VTKDataArrays, name + "IsSet", properties.size(), selection, eType);
		data<double, vtkDoubleArray>(VTKGrid, VTKDataArrays, name + "FixedValue", properties.size(), values, eType);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown element type";
	}
}

void VTK::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	data<double, vtkDoubleArray>(VTKGrid, VTKDataArrays, name, dimension, values, eType);
}

void VTK::mesh(const OutputConfiguration &output, const Mesh &mesh, const std::string &path, ElementType eType)
{
	VTK vtk(output, mesh, path);

	points(vtk.VTKGrid, mesh, [&] (const Point &point, size_t part) { return vtk.shrink(point, part); });
	cells(vtk.VTKGrid, vtk.VTKDataArrays, mesh, eType);
	save(vtk.VTKGrid, path, output);
}

void VTK::fixPoints(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
{
	std::vector<Element*> fixPoints;
	for (size_t p = 0; p < mesh.parts(); p++) {
		fixPoints.insert(fixPoints.end(), mesh.fixPoints(p).begin(), mesh.fixPoints(p).end());
	}

	std::sort(fixPoints.begin(), fixPoints.end());
	Esutils::removeDuplicity(fixPoints);

	VTK vtk(output, mesh, path);

	points(vtk.VTKGrid, vtk.VTKDataArrays, mesh, fixPoints, [&] (const Point &point, size_t part) { return vtk.shrink(point, part); });
	save(vtk.VTKGrid, path, output);
}

void VTK::corners(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
{
	VTK vtk(output, mesh, path);

	points(vtk.VTKGrid, vtk.VTKDataArrays, mesh, mesh.corners(), [&] (const Point &point, size_t part) { return vtk.shrink(point, part); });
	save(vtk.VTKGrid, path, output);
}

void VTK::gluing(const OutputConfiguration &output, const Mesh &mesh, const Constraints &constraints, const std::string &path, size_t dofs)
{
	VTK vtk(output, mesh, "");

	std::vector<Point> cCenter(environment->MPIsize);
	std::vector<Point> sCenters(environment->MPIsize * mesh.parts());

	MPI_Allgather(&vtk._cCenter, sizeof(Point), MPI_BYTE, cCenter.data(), sizeof(Point), MPI_BYTE, MPI_COMM_WORLD);
	MPI_Allgather(vtk._sCenters.data(), mesh.parts() * sizeof(Point), MPI_BYTE, sCenters.data(), mesh.parts() * sizeof(Point), MPI_BYTE, MPI_COMM_WORLD);

	std::vector<size_t> sOffset(mesh.parts());
	std::vector<size_t> eOffset(mesh.parts());

	std::vector<std::vector<std::vector<std::pair<esglobal, Element*> > > > DOF2e(mesh.parts(), std::vector<std::vector<std::pair<esglobal, Element*> > >(dofs));

	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
		eOffset[p] += std::lower_bound(constraints.B1[p].I_row_indices.begin(), constraints.B1[p].I_row_indices.end(), constraints.block[Constraints::DIRICHLET] + 1) - constraints.B1[p].I_row_indices.begin();
	}
	for (size_t n = 0; n < mesh.nodes().size(); n++) {
		for (size_t d = 0; d < mesh.nodes()[n]->domains().size(); d++) {
			size_t p = mesh.nodes()[n]->domains()[d];
			for (size_t dof = 0; dof < dofs; dof++) {
				DOF2e[p][dof].push_back(std::make_pair(mesh.nodes()[n]->DOFIndex(p, dof) + IJVMatrixIndexing, mesh.nodes()[n]));
			}
		}
	}
	cilk_for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t dof = 0; dof < dofs; dof++) {
			std::sort(DOF2e[p][dof].begin(), DOF2e[p][dof].end());
		}
	}

	for (size_t dof = 0; dof < dofs; dof++) {
		std::stringstream ss;
		ss << dof << "dirichlet";
		VTK dirichlet(output, mesh, ss.str());

		std::vector<std::vector<eslocal> > dnodes(mesh.parts());

		cilk_for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = sOffset[p]; i < eOffset[p]; i++) {
				auto it = std::lower_bound(DOF2e[p][dof].begin(), DOF2e[p][dof].end(), constraints.B1[p].J_col_indices[i], [] (const std::pair<esglobal, Element*> &pair, esglobal index) {
					return pair.first < index;
				});
				if (it != DOF2e[p][dof].end() && it->first == constraints.B1[p].J_col_indices[i]) {
					dnodes[p].push_back(it->second->node(0));
				}
			}
		}

		points(dirichlet.VTKGrid, dirichlet.VTKDataArrays, mesh, dnodes, [&] (const Point &point, size_t part) { return vtk.shrink(point, part); });
		save(dirichlet.VTKGrid, ss.str(), output);
	}

	size_t maxLambda = constraints.block[Constraints::DIRICHLET] + constraints.block[Constraints::EQUALITY_CONSTRAINTS] + 1;
	for (size_t p = 0; p < mesh.parts(); p++) {
		sOffset[p] = eOffset[p];
		eOffset[p] = std::lower_bound(constraints.B1[p].I_row_indices.begin() + eOffset[p], constraints.B1[p].I_row_indices.end(), maxLambda) - constraints.B1[p].I_row_indices.begin();
	}

	auto getDomain = [&] (esglobal lambda, size_t exclude) -> size_t {
		for (size_t p = 0; p < mesh.parts(); p++) {
			if (p == exclude) {
				continue;
			}
			auto it = std::lower_bound(constraints.B1subdomainsMap[p].begin(), constraints.B1subdomainsMap[p].end(), lambda - IJVMatrixIndexing);
			if (it != constraints.B1subdomainsMap[p].end() && *it == lambda - IJVMatrixIndexing) {
				return p;
			}
		}
		ESINFO(ERROR) << "Internal error: Broken exporting of gluing matrices";
		return 0;
	};

	for (size_t dof = 0; dof < dofs; dof++) {
		std::stringstream ss;
		ss << dof << "gluing";
		VTK gluing(output, mesh, ss.str());

		std::vector<std::vector<eslocal> > dnodes(mesh.parts());
		std::vector<std::vector<size_t> > indices(mesh.parts());

		std::vector<std::vector<std::vector<std::pair<esglobal, eslocal> > > > CDneighbour(mesh.parts(), std::vector<std::vector<std::pair<esglobal, eslocal> > >(environment->MPIsize));
		std::vector<std::vector<std::pair<eslocal, eslocal> > > CDindex(mesh.parts()); // cluster x domain index
		std::vector<std::vector<std::pair<esglobal, eslocal> > > sBuffer(environment->MPIsize);
		std::vector<std::vector<std::pair<esglobal, eslocal> > > rBuffer(environment->MPIsize);

		cilk_for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = sOffset[p]; i < eOffset[p]; i++) {
				auto it = std::lower_bound(DOF2e[p][dof].begin(), DOF2e[p][dof].end(), constraints.B1[p].J_col_indices[i], [] (const std::pair<esglobal, Element*> &pair, esglobal index) {
					return pair.first < index;
				});
				if (it != DOF2e[p][dof].end() && it->first == constraints.B1[p].J_col_indices[i]) {
					dnodes[p].push_back(it->second->node(0));
					indices[p].push_back(i);

					auto it = std::lower_bound(constraints.B1clustersMap.begin(), constraints.B1clustersMap.end(), constraints.B1[p].I_row_indices[i] - IJVMatrixIndexing, [&] (const std::vector<esglobal> &v, esglobal i) {
						return v[0] < i;
					});
					if (it->size() == 2) { // local gluing
						CDindex[p].push_back(std::make_pair((eslocal)it->back(), (eslocal)getDomain(constraints.B1[p].I_row_indices[i], p)));
					} else { // global gluing
						CDindex[p].push_back(std::make_pair((eslocal)it->back(), -1));
						CDneighbour[p][it->back()].push_back(std::make_pair(constraints.B1[p].I_row_indices[i], p));
					}

				}
			}
		}

		cilk_for (int c = 0; c < environment->MPIsize; c++) {
			for (size_t p = 0; p < mesh.parts(); p++) {
				sBuffer[c].insert(sBuffer[c].end(), CDneighbour[p][c].begin(), CDneighbour[p][c].end());
			}
			std::sort(sBuffer[c].begin(), sBuffer[c].end(), [] (const std::pair<esglobal, eslocal> &p1, const std::pair<esglobal, eslocal> &p2) {
				return p1.first < p2.first;
			});
		}

		std::vector<MPI_Request> req(2 * environment->MPIsize);
		for (int n = 0; n < environment->MPIsize; n++) {
			rBuffer[n].resize(sBuffer[n].size());
			MPI_Isend(sBuffer[n].data(), sizeof(std::pair<esglobal, eslocal>) * sBuffer[n].size(), MPI_BYTE, n, 1, MPI_COMM_WORLD, req.data() + 2 * n);
			MPI_Irecv(rBuffer[n].data(), sizeof(std::pair<esglobal, eslocal>) * rBuffer[n].size(), MPI_BYTE, n, 1, MPI_COMM_WORLD, req.data() + 2 * n + 1);
		}

		MPI_Waitall(2 * environment->MPIsize, req.data(), MPI_STATUSES_IGNORE);

		for (size_t p = 0; p < mesh.parts(); p++) {
			std::vector<int> offsets(environment->MPIsize);
			for (size_t i = 0; i < CDindex[p].size(); i++) {
				if (CDindex[p][i].second == -1) {
					int n = CDindex[p][i].first;
					esglobal lambda = CDneighbour[p][n][offsets[n]++].first;
					auto it = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), lambda, [] (const std::pair<esglobal, eslocal> &p, esglobal l) { return p.first < l; });
					if (it != rBuffer[n].end() && it->first == lambda) {
						CDindex[p][i].second = it->second;
					} else {
						ESINFO(ERROR) << "Internal error: gluing matrices are weird.";
					}
				}
			}
		}

		size_t cSize = 0;
		for (size_t p = 0; p < dnodes.size(); p++) {
			cSize += dnodes[p].size();
		}

		lambdas(gluing.VTKGrid, gluing.VTKDataArrays, mesh, dnodes, [&] (const Point &point, size_t part, size_t index, bool start) -> Point {
			if (start) {
				return vtk.shrink(point, part);
			} else {
				Point p1 = vtk.shrink(point, part);
				Point p2 = vtk.shrink(point, sCenters[CDindex[part][index].first * mesh.parts() + CDindex[part][index].second], cCenter[CDindex[part][index].first]);
				return p1 + (p2 - p1) * constraints.B1duplicity[part][indices[part][index]];
			}
		});

		std::vector<std::vector<double> > values(mesh.parts());
		for (size_t p = 0; p < dnodes.size(); p++) {
			for (size_t n = 0; n < dnodes[p].size(); n++) {
				values[p].push_back(constraints.B1[p].V_values[indices[p][n]]);
			}
		}

		data<double, vtkDoubleArray>(gluing.VTKGrid, gluing.VTKDataArrays, "values", 1, values, ElementType::ELEMENTS);
		gluing.VTKGrid->GetCellData()->SetActiveScalars("values");
		save(gluing.VTKGrid, ss.str(), output);
	}
}

void VTK::finalize()
{
	save(VTKGrid, _path, _output);
}
