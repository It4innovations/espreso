
#include "../vtkxmlbinary.h"

#include <cstring>
#include <fstream>

#include "../../../basis/logging/logging.h"

using namespace espreso::output;

VTKXMLBinary::VTKXMLBinary(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode)
: VTKXML(output, mesh, mode)
{

}

VTKXMLBinary::~VTKXMLBinary()
{

}

struct Base64Decoder {

	static std::string base64char;

	char buffer[3] = { 0, 0, 0 };
	int bufferSize = 0;

	template <typename Ttype>
	void add(const Ttype &data)
	{
		size_t i = 0;
		if (bufferSize) {
			if (sizeof(Ttype) >= 3 - bufferSize) {
				memcpy(buffer + bufferSize, &data, 3 - bufferSize);
				store(buffer, 3);
				i = 3 - bufferSize;
				bufferSize = 0;
			} else {
				memcpy(buffer + bufferSize, &data, sizeof(Ttype));
				bufferSize += sizeof(Ttype);
				i = sizeof(Ttype);
			}
		}

		while (i + 3 < sizeof(Ttype)) {
			store(((char*)&data) + i, 3);
			i += 3;
		}
		if (sizeof(Ttype) % 3 > 0) {
			memset(buffer, 0, 3);
			memcpy(buffer, ((char*)&data) + i, sizeof(Ttype) - i);
			bufferSize = sizeof(Ttype) - i;
		}
	}

	void finalize()
	{
		if (bufferSize) {
			store(buffer, bufferSize);
		}
		bufferSize = 0;
	}

	Base64Decoder(std::ofstream &os): _os(os) {};

	~Base64Decoder()
	{
		finalize();
	}

private:
	std::ofstream &_os;

	void store(char *data, int len)
	{
		_os << base64char[((data[0] >> 2) & 0x3F)];
		_os << base64char[(((data[0] << 4) & 0x30) | ((data[1] >> 4) & 0x0F))];
		_os << ((len >= 2 ? base64char[((data[1] << 2) & 0x3C) | ((data[2] >> 6) & 0x03)] : '='));
		_os << ((len == 3 ? base64char[(data[2] & 0x3F)] : '='));
	}
};

std::string Base64Decoder::base64char = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";


template <typename TType>
static void storeData(std::ofstream &os, const std::string &type, const std::string &name, size_t components, const std::vector<TType> &data)
{
	os << "    <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"binary\" NumberOfComponents=\"" << components << "\">\n";
	os << "      ";

	int size = data.size();
	Base64Decoder decoder(os);
	decoder.add((int)(size * sizeof(TType)));
	decoder.finalize();
	for (size_t i = 0; i < data.size(); i++) {
		decoder.add(data[i]);
	}
	decoder.finalize();

	os << "\n";
	os << "    </DataArray>\n";
}

void VTKXMLBinary::storePointData(const std::string &name, size_t components, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, components, data);
}

void VTKXMLBinary::storePointData(const std::string &name, size_t components, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, components, data);
}

void VTKXMLBinary::storePointData(const std::string &name, size_t components, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, components, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t components, const std::vector<int> &data)
{
	storeData(*_os, "Int32", name, components, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t components, const std::vector<long> &data)
{
	storeData(*_os, "Int64", name, components, data);
}

void VTKXMLBinary::storeCellData(const std::string &name, size_t components, const std::vector<double> &data)
{
	storeData(*_os, "Float64", name, components, data);
}



