
#include <fstream>

#include "../vtkxmlbinary.h"

using namespace espreso::output;

VTKXMLBinary::VTKXMLBinary(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: VTKXML(output, mesh, path)
{

}

void VTKXMLBinary::store(std::ostream &os, const std::vector<eslocal> &data)
{
	os.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(eslocal));
}

void VTKXMLBinary::store(std::ostream &os, const std::vector<double> &data)
{
	os.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
}



