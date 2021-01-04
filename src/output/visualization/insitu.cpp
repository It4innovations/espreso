
#include "insitu.h"

#include "esinfo/ecfinfo.h"

#include <unistd.h>
#include "wrappers/catalyst/w.catalyst.h"

using namespace espreso;

InSitu::InSitu()
: _catalyst(NULL)
{

}

InSitu::~InSitu()
{
	if (_catalyst != NULL) {
		delete _catalyst;
	}
}

void InSitu::updateMesh()
{
	_catalyst = new Catalyst();
}

void InSitu::updateSolution()
{
	_catalyst->update();
	sleep(info::ecf->output.catalyst_sleep_time);
}




