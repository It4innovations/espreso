
#include "../vtkxmlascii.h"

#include <fstream>


using namespace espreso::output;

VTKXMLASCII::VTKXMLASCII(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: VTKXML(output, mesh, path)
{

}

VTKXMLASCII::~VTKXMLASCII()
{

}

template <typename TType>
static void storeData(std::ofstream &os, const std::string &type, const std::string &name, size_t elements, const std::vector<std::vector<TType> > &data)
{
	size_t size = 0;
	for (size_t i = 0; i < data.size(); i++) {
		size += data[i].size();
	}
	os << "    <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\" NumberOfComponents=\"" << size / elements << "\">\n";
	os << "      ";
	for (size_t i = 0; i < data.size(); i++) {
		for (size_t j = 0; j < data[i].size(); j++) {
			os << data[i][j] << " ";
		}
	}
	os << "\n";
	os << "    </DataArray>\n";
}

template <typename TType>
static void storeData(std::ofstream &os, const std::string &type, const std::string &name, size_t elements, const std::vector<TType> &data)
{
	os << "    <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\" NumberOfComponents=\"" << data.size() / elements << "\">\n";
	os << "      ";
	for (size_t i = 0; i < data.size(); i++) {
		os << data[i] << " ";
	}
	os << "\n";
	os << "    </DataArray>\n";
}

void VTKXMLASCII::storePointData(const std::string &name, size_t points, const std::vector<std::vector<int> > &data)
{
	storeData(*_os, "Int32", name, points, data);
}

void VTKXMLASCII::storePointData(const std::string &name, size_t points, const std::vector<std::vector<long> > &data)
{
	storeData(*_os, "Int64", name, points, data);
}

void VTKXMLASCII::storePointData(const std::string &name, size_t points, const std::vector<std::vector<double> > &data)
{

	storeData(*_os, "Float64", name, points, data);
}

void VTKXMLASCII::storePointData(const std::string &name, size_t points, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, points, data);
}

void VTKXMLASCII::storePointData(const std::string &name, size_t points, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, points, data);
}

void VTKXMLASCII::storePointData(const std::string &name, size_t points, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, points, data);
}


void VTKXMLASCII::storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<int> > &data)
{
	storeData(*_os, "Int32", name, cells, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<long> > &data)
{
	storeData(*_os, "Int64", name, cells, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<double> > &data)
{
	storeData(*_os, "Float64", name, cells, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t cells, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, cells, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t cells, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, cells, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t cells, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, cells, data);
}



