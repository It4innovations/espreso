
#include "instance.h"

Instance::Instance(int rank, int size)
	: _mesh(rank, size), _rank(rank), _size(size)
{

}



