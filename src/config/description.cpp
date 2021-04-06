
#include "description.h"
#include "configuration.h"

using namespace espreso;

ECFDescription::ECFDescription()
{
	ecfdescription = new ECFObject();
}

ECFDescription::~ECFDescription()
{
	delete ecfdescription;
}
