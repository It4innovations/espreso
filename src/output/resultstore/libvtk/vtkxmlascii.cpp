
#include "../vtkxmlascii.h"

#include <fstream>

#include "vtkXMLUnstructuredGridWriter.h"

using namespace espreso::output;

VTKXMLASCII::VTKXMLASCII(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: VTKXML(output, mesh, path)
{

}

void VTKXMLASCII::store(std::ostream &os, const std::vector<eslocal> &data)
{

}

void VTKXMLASCII::store(std::ostream &os, const std::vector<double> &data)
{

}



