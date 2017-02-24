
// Dummy Catalyst file
// Do nothing

#include "../catalyst.h"
#include "../../../basis/logging/logging.h"

using namespace espreso::store;

Catalyst::Catalyst(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
: VTK(output, mesh, path), processor(NULL), dataDescription(NULL)
{
	ESINFO(ALWAYS) << Info::TextColor::YELLOW << "ESPRESO not supports Catalyst - re-compile ESPRESO with VTK library.";
}

Catalyst::~Catalyst()
{

}

void Catalyst::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType) { }
void Catalyst::finalize() { }
