
#include "parser.h"

#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"

using namespace espreso;

OpenFOAMParser::OpenFOAMParser(const char* begin, const char* end)
: begin(begin), end(end)
{

}

OpenFOAMSeparateParser::OpenFOAMSeparateParser(const char *begin, const char *end)
: OpenFOAMParser(begin, end)
{
	while (this->begin != this->end && *this->begin++ != '(');
	while (this->end != this->begin && *(--this->end) != ')');
}

OpenFOAMCollectiveParser::OpenFOAMCollectiveParser(const char *begin, const char *end)
: OpenFOAMParser(begin, end)
{
	int found, min, max;
	const char *p;

	p = begin;
	while (*p != '(' && p < end) { ++p; };
	if (p == end) {
		found = info::mpi::size;
	} else {
		found = info::mpi::rank;
	}

	Communication::allReduce(&found, &min, 1, MPI_INT, MPI_MIN);
	if (info::mpi::rank == min) {
		this->begin = p + 1;
	}
	if (info::mpi::rank < min) {
		this->begin = p;
	}

	p = end;
	while (*p != ')' && begin < p) { --p; }
	if (p == begin) {
		found = 0;
	} else {
		found = info::mpi::rank;
	}

	Communication::allReduce(&found, &max, 1, MPI_INT, MPI_MAX);
	if (max <= info::mpi::rank) {
		this->end = p;
	}
}





