
#include "../catalyst.h"

#include "../../../basis/logging/logging.h"

using namespace espreso;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode)
: VTKXMLASCII(output, mesh, mode & ~MeshInfo::PREPARE), _processor(NULL), _dataDescription(NULL), _fieldData(NULL)
{
	ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Link VTK library to support Paraview Catalyst.";
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

void Catalyst::storeFETIData(const Step &step, const Instance &instance)
{

}

void Catalyst::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{

}



