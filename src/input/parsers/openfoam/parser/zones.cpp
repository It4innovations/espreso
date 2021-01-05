
#include "zones.h"

#include "input/parsers/openfoam/openfoam.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"

using namespace espreso;

OpenFOAMZones::OpenFOAMZones(InputFile &pfile)
: OpenFOAMCollectiveParser(pfile.begin, pfile.end), _pfile(pfile)
{

}

int OpenFOAMZones::getZones()
{
	int zones = 0, root = 0, rank = 0;
	const char *c = begin;
	if (*(c - 1) == '(') {
		c -= 3;
		while (*c != '\n') { c--; } // go before number of boundaries
		rank = info::mpi::rank;
		zones = readInteger(c);
	}

	Communication::allReduce(&rank, &root, 1, MPI_INT, MPI_SUM);
	Communication::broadcast(&zones, 1, MPI_INT, root);
	return zones;
}

void OpenFOAMZones::synchronize(int zones, std::vector<char> &names, std::vector<size_t> &offsets)
{
	names.resize(zones * 80);

	std::vector<char> mynames(zones * 80);
	memset(mynames.data(), '\0', mynames.size());
	int nbrackets = 0;
	int noffsets = 0;
	std::vector<const char*> brackets;

	if (begin != end) {
		const char *c = begin + 1, *line = begin;
		while (c != end) {
			if (*c == '\n') {
				line = c + 1;
			}
			if (*c == '{') {
				if (!StringCompare::caseInsensitiveEq(readString(line), "flipMap")) {
					brackets.push_back(c);
				}
			}
			if (*c == '(' || *c == ')') {
				offsets.push_back(_pfile.distribution[info::mpi::rank] + (c - _pfile.begin));
			}
			++c;
		}
		if (*c == '{') { // if '{' is the first of the next process, I have to parse the name
			brackets.push_back(c);
		}
	}

	nbrackets = brackets.size();
	noffsets = offsets.size();
	Communication::exscan(nbrackets);
	Communication::exscan(noffsets);
	std::vector<size_t> _offset(noffsets);
	_offset.insert(_offset.end(), offsets.begin(), offsets.end());
	_offset.resize(2 * zones);
	offsets.resize(2 * zones);
	Communication::allReduce(_offset.data(), offsets.data(), _offset.size(), MPITools::getType<size_t>().mpitype, MPI_SUM);

	for (size_t i = 0; i < brackets.size(); i++) {
		const char *name = brackets[i];
		while (*name-- != '\n');
		while (*name != '\n') { name--; }
		++name;
		std::string zonename = readString(name);
		memcpy(mynames.data() + 80 * (noffsets + i), zonename.data(), zonename.size());
	}

	Communication::allReduce(mynames.data(), names.data(), names.size(), MPI_CHAR, MPI_SUM);
}

void OpenFOAMZones::readData(std::vector<esint> &indices, size_t begin, size_t end)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	if (begin + 1 > _pfile.distribution[info::mpi::rank + 1]) {
		return;
	}
	if (end - 1 < _pfile.distribution[info::mpi::rank]) {
		return;
	}

	begin = std::max(begin + 1, _pfile.distribution[info::mpi::rank]);
	end = std::min(end - 1, _pfile.distribution[info::mpi::rank + 1]);
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);

	std::vector<std::vector<esint> > data(threads);
	size_t offset = begin - _pfile.distribution[info::mpi::rank];

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;

		const char *c = _pfile.begin + offset + tdistribution[t];
		if (c > _pfile.begin) {
			while (c < _pfile.end && *(c - 1) != '\n') { ++c; }
		}
		while (c < _pfile.begin + offset + tdistribution[t + 1]) {
			tdata.push_back(readInteger(c));
			c += 1; // skip '\n'
		}

		data[t].swap(tdata);
	}

	for (size_t t = 0; t < threads; t++) {
		indices.insert(indices.end(), data[t].begin(), data[t].end());
	}
	std::sort(indices.begin(), indices.end());
}

bool OpenFOAMZones::readPoints(OpenFOAMData &data)
{
	int zones = getZones();
	if (zones == 0) {
		return true;
	}

	std::vector<char> names;
	std::vector<size_t> offsets;
	synchronize(zones, names, offsets);

	for (int i = 0; i < zones; i++) {
		std::string name(names.data() + 80 * i);
		auto &indices = data.nregions[name];
		readData(indices, offsets[2 * i], offsets[2 * i + 1]);
	}

	return true;
}

bool OpenFOAMZones::readFaces(OpenFOAMData &data)
{
	int zones = getZones();
	if (zones == 0) {
		return true;
	}

	std::vector<char> names;
	std::vector<size_t> offsets;
	synchronize(zones, names, offsets);

	for (int i = 0; i < zones; i++) {
		std::string name(names.data() + 80 * i);
		auto &indices = data.eregions[name];
		readData(indices, offsets[2 * i], offsets[2 * i + 1]);
	}

	return true;
}

bool OpenFOAMZones::readCells(OpenFOAMData &data)
{
	int zones = getZones();
	if (zones == 0) {
		return true;
	}

	std::vector<char> names;
	std::vector<size_t> offsets;
	synchronize(zones, names, offsets);

	for (int i = 0; i < zones; i++) {
		std::string name(names.data() + 80 * i);
		auto &indices = data.eregions[OpenFOAMLoader::cellprefix + name];
		readData(indices, offsets[2 * i], offsets[2 * i + 1]);
	}

	return true;
}






