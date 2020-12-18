
#include "nblock.h"

#include "basis/containers/point.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"

#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"

#include "input/parsers/ansyscdb/ansyscdb.h"

using namespace espreso;

NBlock::NBlock()
: NUMFIELD(-1), Solkey(0), NDMAX(-1), NDSEL(-1),
  lineSize(-1), lineEndSize(-1),
  indexSize(-1), indexLength(-1), valueSize(-1), valueLength(-1)
{

}

NBlock& NBlock::parse(const char* begin)
{
	auto error = [&] (std::string &line) {
		eslog::error("Workbench parse error: unknown format of NBLOCK: '%s'\n", line.c_str());
	};

	std::string commandLine = Parser::getLine(begin);
	std::string descriptionLine = Parser::getLine(begin + commandLine.size());

	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	std::vector<std::string> description = Parser::split(descriptionLine.substr(1, descriptionLine.size() - 2 - lineEndSize), ",");
	if (description.size() != 2) {
		error(descriptionLine);
	}
	std::vector<std::string> indices = Parser::split(description[0], "i");
	std::vector<std::string> values = Parser::split(description[1], "e");

	if (indices.size() != 2 || values.size() < 2) {
		error(descriptionLine);
	}

	switch (command.size()) {
	case 5:
		NDSEL = command[4].size() ? std::stol(command[4]) : -1;
	case 4:
		NDMAX = command[3].size() ? std::stol(command[3]) : -1;
	case 3:
		Solkey = command[2].size();
	case 2:
		NUMFIELD = std::stoi(command[1]);
		break;
	default:
		error(commandLine);
	}

	indexSize = std::stoi(indices[0]);
	indexLength = std::stoi(indices[1]);
	valueSize = std::stoi(values[0]);
	valueLength = std::stoi(values[1]);

	lineSize = indexSize * indexLength + valueSize * valueLength + lineEndSize; // just guess since exceptions allowed by Ansys

	WorkbenchParser::fillIndices(begin, begin + commandLine.size() + descriptionLine.size());
	return *this;
}

bool NBlock::readData(AnsysCDBData &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, last - first);

	std::vector<std::vector<esint> > tids(threads);
	std::vector<std::vector<Point> > tcoords(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<char> index(indexLength + 1);
		std::vector<esint> ids;
		std::vector<Point> coords;
		ids.reserve(1.1 * (tdistribution[t + 1] - tdistribution[t]) / lineSize);
		coords.reserve(1.1 * (tdistribution[t + 1] - tdistribution[t]) / lineSize);

		auto data = first + tdistribution[t];
		if (tdistribution[t] && *(data - 1) != '\n') {
			while (data < first + tdistribution[t + 1] && *data++ != '\n'); // start at new line
		}

		auto getindex = [&] () {
			memcpy(index.data(), data, indexLength);
			data += indexLength;
			return atol(index.data()) - 1;
		};
		auto getvalue = [&] () {
			if (*data != '\r' && *data != '\n') {
				memcpy(value.data(), data, valueLength);
				data += valueLength;
				return atof(value.data());
			}
			return .0;
		};

		while (data < first + tdistribution[t + 1]) {
			ids.push_back(getindex());
			data += (indexSize - 1) * indexLength; // skip solid and line indices
			double x = getvalue();
			double y = getvalue();
			double z = getvalue();
			coords.push_back(Point(x, y, z));
			while (data < first + tdistribution[t + 1] && *data++ != '\n'); // finish line (skip rotation)
		}

		tcoords[t].swap(coords);
		tids[t].swap(ids);
	}

	size_t totalsize = 0;
	for (size_t t = 0; t < threads; t++) {
		totalsize += tids[t].size();
	}
	mesh.nIDs.reserve(mesh.nIDs.size() + totalsize);
	mesh.coordinates.reserve(mesh.coordinates.size() + totalsize);

	for (size_t t = 0; t < threads; t++) {
		mesh.nIDs.insert(mesh.nIDs.end(), tids[t].begin(), tids[t].end());
		mesh.coordinates.insert(mesh.coordinates.end(), tcoords[t].begin(), tcoords[t].end());
	}
	return true;
}
