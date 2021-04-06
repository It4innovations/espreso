
#include "parser.h"
#include "blockend.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"

#include <iostream>
#include <algorithm>

using namespace espreso;

size_t AbaqusParser::offset = 0;
const char* AbaqusParser::begin = NULL;
const char* AbaqusParser::end = NULL;

void AbaqusParser::fillIndices(const char* header, const char* data)
{
	this->header = offset + header - begin;
	this->first = offset + data - begin;
}

void AbaqusParser::fillIndices(const char* header, const char* first, const char* last)
{
	fillIndices(header, first);
	this->last = offset + last - begin;
}

const char* AbaqusParser::getFirst() const
{
	if (fRank <= info::mpi::rank && info::mpi::rank <= lRank) {
		if (fRank == info::mpi::rank) {
			return begin + first - offset;
		} else {
			return begin;
		}
	}
	return begin;
}

const char* AbaqusParser::getLast() const
{
	if (fRank <= info::mpi::rank && info::mpi::rank <= lRank) {
		if (lRank == info::mpi::rank) {
			return begin + last - offset;
		} else {
			return end;
		}
	}
	return begin;
}

void AbaqusParser::fillDistribution(std::vector<BlockFinish> &blocksFinishs, std::vector<size_t> &distribution)
{
	if (last == (size_t)-1) {
		last = std::lower_bound(blocksFinishs.begin(), blocksFinishs.end(), first, [] (BlockFinish &b, size_t first) { return b.first < first; })->first;
	}

	fRank = std::lower_bound(distribution.begin(), distribution.end(), first + 1) - distribution.begin() - 1;
	lRank = std::lower_bound(distribution.begin(), distribution.end(), last + 1) - distribution.begin() - 1;
}

std::string AbaqusParser::command() const
{
	std::string cmd;
	const char *p = getFirst();
	while (*p != '\n') {
		cmd += *p++;
	}
	return cmd;
}

void AbaqusParser::print(const char* data)
{
	std::cout << first << "(" << fRank << ") -> " << last << "(" << lRank << ")\n";
	if (fRank == info::mpi::rank) {
		const char *p = data + first - offset;
		while (*p != '\n') {
			std::cout << *p++;
		}
		std::cout << "\n";
	}
	if (lRank == info::mpi::rank) {
		const char *p = data + last - offset - 1;
		while (*(p - 1) != '\n') { --p; }
		while (*p != '\n') {
			std::cout << *p++;
		}
		std::cout << "\n";
	}
}


