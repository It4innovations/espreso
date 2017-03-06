
#include "../vtkxmlbinary.h"

#include <fstream>

#include "../../../basis/logging/logging.h"

using namespace espreso::output;

VTKXMLBinary::VTKXMLBinary(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: VTKXML(output, mesh, path)
{
	ESINFO(GLOBAL_ERROR) << "VTK XML binaty format is not working without linked VTK library.";
}

VTKXMLBinary::~VTKXMLBinary()
{

}

template <typename TType>
static void storeData(std::ofstream &os, const std::string &type, const std::string &name, size_t elements, const std::vector<std::vector<TType> > &data)
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		size += data[i].size();
	}
	os << "    <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"binary\" NumberOfComponents=\"" << elements / size << "\">\n";
	os << "      ";
	for (size_t i = 0; i < data.size(); i++) {
		os.write(reinterpret_cast<const char*>(data[i].data()), data[i].size() * sizeof(TType));
	}
	os << "\n";
	os << "    </DataArray>\n";
}

template <typename TType>
static void storeData(std::ofstream &os, const std::string &type, const std::string &name, size_t elements, const std::vector<TType> &data)
{
	os << "    <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"binary\" NumberOfComponents=\"" << elements / data.size() << "\">\n";
	os << "      ";
	for (size_t i = 0; i < data.size(); i++) {
		os.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(TType));
	}
	os << "\n";
	os << "    </DataArray>\n";
}

void VTKXMLBinary::storePointData(const std::string &name, size_t points, const std::vector<std::vector<int> > &data)
{
	storeData(*_os, "Int32", name, points, data);
}

void VTKXMLBinary::storePointData(const std::string &name, size_t points, const std::vector<std::vector<long> > &data)
{
	storeData(*_os, "Int64", name, points, data);
}

void VTKXMLBinary::storePointData(const std::string &name, size_t points, const std::vector<std::vector<double> > &data)
{

	storeData(*_os, "Float64", name, points, data);
}

void VTKXMLBinary::storePointData(const std::string &name, size_t points, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, points, data);
}

void VTKXMLBinary::storePointData(const std::string &name, size_t points, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, points, data);
}

void VTKXMLBinary::storePointData(const std::string &name, size_t points, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, points, data);
}


void VTKXMLBinary::storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<int> > &data)
{
	storeData(*_os, "Int32", name, cells, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<long> > &data)
{
	storeData(*_os, "Int64", name, cells, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<double> > &data)
{
	storeData(*_os, "Float64", name, cells, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t cells, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, cells, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t cells, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, cells, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t cells, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, cells, data);
}



