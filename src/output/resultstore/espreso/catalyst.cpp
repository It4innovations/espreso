
#include "../catalyst.h"

#include "../../../basis/logging/logging.h"

using namespace espreso::output;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: VTKXMLASCII(output, mesh, path), _processor(NULL), _dataDescription(NULL), _fieldData(NULL)
{
	ESINFO(GLOBAL_ERROR) << "Link VTK library to support Paraview Catalyst.";
}

Catalyst::~Catalyst()
{

}

void Catalyst::storeSettings(const Step &step)
{

}

void Catalyst::storeSettings(size_t steps)
{

}

void Catalyst::storeSettings(const std::vector<size_t> &steps)
{

}

void Catalyst::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{

}

void Catalyst::finalize()
{

}



