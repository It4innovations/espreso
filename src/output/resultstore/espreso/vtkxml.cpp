
#include <fstream>
#include <algorithm>

#include "../../../configuration/environment.h"
#include "../../../assembler/solution.h"
#include "../../../mesh/structures/elementtypes.h"

#include "../vtkxmlascii.h"

using namespace espreso::output;

VTKXML::VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path)
{
	preprocessing();
	size_t offset = 0;
	std::for_each(_elementsNodes.begin(), _elementsNodes.end(), [&] (eslocal &nodes) { nodes = offset += nodes; });
}

void VTKXML::regionPreprocessing(const espreso::Region &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	ResultStore::regionPreprocessing(region, coordinates, elementsTypes, elementsNodes, elements);
	size_t offset = 0;
	std::for_each(elementsNodes.begin(), elementsNodes.end(), [&] (eslocal &nodes) { nodes = offset += nodes; });
}

void VTKXML::storeMesh(std::ofstream &os, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
	os << "  <Points>\n";
	os << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << format() << "\">\n";
	os << "      "; store(os, coordinates); os << "\n";
	os << "    </DataArray>\n";
	os << "  </Points>\n";

	os << "  <Cells>\n";
	os << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"" << format() << "\">\n";
	os << "      "; store(os, elements); os << "\n";
	os << "    </DataArray>\n";
	os << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"" << format() << "\">\n";
	os << "      "; store(os, elementsNodes); os << "\n";
	os << "    </DataArray>\n";
	os << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"" << format() << "\">\n";
	os << "      "; store(os, elementsTypes); os << "\n";
	os << "    </DataArray>\n";
	os << "  </Cells>\n";
}

static void xmlHeat(std::ofstream &os, size_t points, size_t cells)
{
	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "<UnstructuredGrid>\n";
	os << "<Piece NumberOfPoints=\"" << points << "\" NumberOfCells=\"" << cells << "\">\n";
}

static void xmlTail(std::ofstream &os)
{
	os << "</Piece>\n";
	os << "</UnstructuredGrid>\n";
	os << "</VTKFile>\n";
}

void VTKXML::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data)
{
	std::ofstream os;
	os.open((name + ".vtu").c_str(), std::ios::out | std::ios::trunc);

	xmlHeat(os, coordinates.size() / 3, elementsTypes.size());
	storeMesh(os, coordinates, elementsTypes, elementsNodes, elements);

	os << "  <PointData>\n";
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		os << "    <DataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() / (coordinates.size() / 3) << "\">\n";
		os << "      "; store(os, *it->second); os << "\n";
		os << "    </DataArray>\n";
	}

	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		os << "    <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() / (coordinates.size() / 3) << "\">\n";
		os << "      "; store(os, *it->second); os << "\n";
		os << "    </DataArray>\n";
	}
	os << "  </PointData>\n";

	os << "  <CellData>\n";
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		os << "    <DataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() / elementsTypes.size() << "\">\n";
		os << "      "; store(os, *it->second); os << "\n";
		os << "    </DataArray>\n";
	}

	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		os << "    <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second->size() / elementsTypes.size() << "\">\n";
		os << "      "; store(os, *it->second); os << "\n";
		os << "    </DataArray>\n";
	}
	os << "  </CellData>\n";

	xmlTail(os);
	os.close();
}

void VTKXML::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution)
{
	std::ofstream os;
	os.open((name + ".vtu").c_str(), std::ios::out | std::ios::trunc);

	xmlHeat(os, coordinates.size() / 3, elementsTypes.size());
	storeMesh(os, coordinates, elementsTypes, elementsNodes, elements);

	os << "  <PointData>\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::NODES) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= (coordinates.size() / 3);
			os << "    <DataArray type=\"Float64\" Name=\"" << solution[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << components  << "\">\n";
			os << "      ";
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { store(os, part); });;
			os << "\n";
			os << "    </DataArray>\n";
		}
	}
	os << "  </PointData>\n";

	os << "  <CellData>\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::ELEMENTS) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= elementsTypes.size();
			os << "    <DataArray type=\"Float64\" Name=\"" << solution[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << components << "\">\n";
			os << "      ";
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { store(os, part); });;
			os << "\n";
			os << "    </DataArray>\n";
		}
	}
	os << "  </CellData>\n";

	xmlTail(os);
	os.close();
}

void VTKXML::savePVTU(const std::string &root, const std::string &name, const std::vector<std::pair<std::string, size_t> > &pointDataInt, const std::vector<std::pair<std::string, size_t>> &pointDataDouble, const std::vector<std::pair<std::string, size_t>> &cellDataInt, const std::vector<std::pair<std::string, size_t>> &cellDataDouble)
{
	std::ofstream os;

	os.open(root + name + ".pvtu", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
	os << "\n";
	os << "    <PPoints>\n";
	os << "      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << format() << "\"/>\n";
	os << "    </PPoints>\n";
	os << "\n";
	os << "    <PPointData>\n";
	for (size_t i = 0; i < pointDataInt.size(); i++) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << pointDataInt[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << pointDataInt[i].second << "\"/>\n";
	}
	for (size_t i = 0; i < pointDataDouble.size(); i++) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << pointDataDouble[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << pointDataDouble[i].second << "\"/>\n";
	}
	os << "    </PPointData>\n";
	os << "\n";
	os << "    <PCellData>\n";
	for (size_t i = 0; i < cellDataInt.size(); i++) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << cellDataInt[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << cellDataInt[i].second << "\"/>\n";
	}
	for (size_t i = 0; i < cellDataDouble.size(); i++) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << cellDataDouble[i].first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << cellDataDouble[i].second << "\"/>\n";
	}
	os << "    </PCellData>\n";
	os << "\n";
	for (int r = 0; r < environment->MPIsize; r++) {
		os << "    <Piece Source=\"" << r << "/" << name << ".vtu\"/>\n";
	}
	os << "  </PUnstructuredGrid>\n";
	os << "</VTKFile>\n";
}

void VTKXML::composeClusters(const std::string &root, const std::string &name, const DataArrays &data)
{
	std::vector<std::pair<std::string, size_t> > pointDataInt, pointDataDouble, cellDataInt, cellDataDouble;

	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		cellDataInt.push_back(std::make_pair(it->first, it->second->size() / _elementsTypes.size()));
	}
	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		cellDataDouble.push_back(std::make_pair(it->first, it->second->size() / _elementsTypes.size()));
	}
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		pointDataInt.push_back(std::make_pair(it->first, it->second->size() / (_coordinates.size() / 3)));
	}
	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		pointDataDouble.push_back(std::make_pair(it->first, it->second->size() / (_coordinates.size() / 3)));
	}
	savePVTU(root, name, pointDataInt, pointDataDouble, cellDataInt, cellDataDouble);
}

void VTKXML::composeClusters(const std::string &root, const std::string &name, const std::vector<Solution*> &solution)
{
	std::vector<std::pair<std::string, size_t> > pointData, cellData;

	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::ELEMENTS) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= _elementsTypes.size();
			cellData.push_back(std::make_pair(solution[i]->name, components));
		}
		if (solution[i]->eType == ElementType::NODES) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= (_coordinates.size() / 3);
			pointData.push_back(std::make_pair(solution[i]->name, components));
		}
	}
	savePVTU(root, name, {}, pointData, {}, cellData);
}

void VTKXML::composeRegions(const std::string &name, const std::vector<std::string> &names)
{
	std::ofstream os;

	os.open(name + ".vtm", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "<vtkMultiBlockDataSet>\n";
	for (size_t i = 0; i < names.size(); i++) {
		os << "  <DataSet index=\"" << i << "\" name=\"" << names[i] << "\" file=\"" << names[i] << ".pvtu\"> </DataSet>\n";
	}
	os << "</vtkMultiBlockDataSet>\n";
	os << "</VTKFile>\n";

}



