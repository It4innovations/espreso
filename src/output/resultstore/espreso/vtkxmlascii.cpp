
#include "../vtkxmlascii.h"

#include <fstream>


using namespace espreso::output;

VTKXMLASCII::VTKXMLASCII(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: VTKXML(output, mesh, path)
{

}

void VTKXMLASCII::store(std::ostream &os, const std::vector<eslocal> &data)
{
	for (size_t i = 0; i < data.size(); i++) {
		os << data[i] << " ";
	}
}

void VTKXMLASCII::store(std::ostream &os, const std::vector<double> &data)
{
	for (size_t i = 0; i < data.size(); i++) {
		os << data[i] << " ";
	}
}



