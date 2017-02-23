
// Dummy Catalyst file
// Do nothing

#include "../catalyst.h"
#include "../../../basis/logging/logging.h"

using namespace espreso::store;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: ResultStore(output, &mesh, path), processor(NULL), VTKGrid(NULL), dataDescription(NULL)
{
	ESINFO(ALWAYS) << Info::TextColor::YELLOW << "ESPRESO not supports Catalyst - re-compile ESPRESO with VTK library.";
}

Catalyst::~Catalyst()
{

}

void Catalyst::storeGeometry(size_t timeStep) { }
void Catalyst::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType) { }
void Catalyst::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType) { }
void Catalyst::finalize() { }
