
#include "instance.h"

Instance::Instance(int rank, int size)
	: _rank(rank), _size(size),_localBoundaries(_mesh), _globalBoundaries(_mesh)
{

}



