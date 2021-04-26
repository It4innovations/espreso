
#include "bodystore.h"

#include "store.h"
#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"


using namespace espreso;

BodyStore::BodyStore()
{
	size = 1;
	totalSize = info::mpi::size;
}

BodyStore::~BodyStore()
{

}
