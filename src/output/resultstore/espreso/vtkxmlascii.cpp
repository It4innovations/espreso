
#include "../vtkxmlascii.h"

#include <fstream>


using namespace espreso;

VTKXMLASCII::VTKXMLASCII(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode)
: VTKXML(output, mesh, mode)
{

}

VTKXMLASCII::~VTKXMLASCII()
{

}

template <typename TType>
static void storeData(std::ofstream &os, const std::string &type, const std::string &name, size_t components, const std::vector<TType> &data)
{
	os << "    <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\" NumberOfComponents=\"" << components << "\">\n";
	os << "      ";
	for (size_t i = 0; i < data.size(); i++) {
		os << data[i] << " ";
	}
	os << "\n";
	os << "    </DataArray>\n";
}

void VTKXMLASCII::storePointData(const std::string &name, size_t components, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, components, data);
}

void VTKXMLASCII::storePointData(const std::string &name, size_t components, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, components, data);
}

void VTKXMLASCII::storePointData(const std::string &name, size_t components, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, components, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t components, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, components, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t components, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, components, data);
}

void VTKXMLASCII::storeCellData(const std::string &name, size_t components, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, components, data);
}



