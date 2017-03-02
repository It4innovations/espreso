
#include <fstream>
#include <algorithm>

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

void VTKXML::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data)
{
	std::ofstream os;

	os.open((name + ".vtu").c_str(), std::ios::out | std::ios::trunc);

	size_t coordinateSize = coordinates.size() / 3;

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "<UnstructuredGrid>\n";
	os << "<Piece NumberOfPoints=\"" << coordinates.size() / 3 << "\" NumberOfCells=\"" << elementsTypes.size() << "\">\n";

	os << "  <Points>\n";
	os << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << format() << "\">\n";
	os << "      "; store(os, coordinates);
	os << "    </DataArray>\n";
	os << "  </Points>\n";

	os << "  <Cells>\n";
	os << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"" << format() << "\">\n";
	os << "      "; store(os, elements);
	os << "    </DataArray>\n";
	os << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"" << format() << "\">\n";
	os << "      "; store(os, elementsNodes);
	os << "    </DataArray>\n";
	os << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"" << format() << "\">\n";
	os << "      "; store(os, elementsTypes);
	os << "    </DataArray>\n";
	os << "  </Cells>\n";

	os << "  <PointData>\n";
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		os << "    <DataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\">\n";
		os << "      "; store(os, *it->second);
		os << "    </DataArray>\n";
	}

	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		os << "    <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\">\n";
		os << "      "; store(os, *it->second);
		os << "    </DataArray>\n";
	}
	os << "  </PointData>\n";

	os << "  <CellData>\n";
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		os << "    <DataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\">\n";
		os << "      "; store(os, *it->second);
		os << "    </DataArray>\n";
	}

	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		os << "    <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\">\n";
		os << "      "; store(os, *it->second);
		os << "    </DataArray>\n";
	}
	os << "  </CellData>\n";


	os << "</Piece>";
	os << "</UnstructuredGrid>";
	os << "</VTKFile>";
}

void VTKXML::compose(const std::string &name, const std::vector<std::string> &names)
{
	std::ofstream os;

	os.open(name + ".vtm", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "<vtkMultiBlockDataSet>\n";
	for (size_t i = 0; i < names.size(); i++) {
		os << "  <DataSet index=\"" << i << "\" name=\"" << names[i] << "\" file=\"" << names[i] << ".vtu\"> </DataSet>\n";
	}
	os << "</vtkMultiBlockDataSet>\n";
	os << "</VTKFile>\n";

}



