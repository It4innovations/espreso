
#include "cmblock.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"

#include <algorithm>

using namespace espreso;

CMBlock::CMBlock()
: entity(Entity::UNKNOWN), NUMITEMS(-1),
  lineSize(-1), lineEndSize(-1),
  valueSize(-1), valueLength(-1)
{
	memset(name, '\0', MAX_NAME_SIZE);
}

CMBlock& CMBlock::parse(const char* begin)
{
	auto error = [&] (std::string &line) {
		eslog::error("Workbench parse error: unknown format of CMBLOCK: %s\n", line.c_str());
	};

	std::string commandLine = Parser::getLine(begin);
	std::string descriptionLine = Parser::getLine(begin + commandLine.size());

	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	std::vector<std::string> description = Parser::split(descriptionLine.substr(1, descriptionLine.size() - 2 - lineEndSize), "i");
	if (command.size() != 4) {
		error(commandLine);
	}
	if (description.size() != 2) {
		error(descriptionLine);
	}

	command[1] = Parser::strip(command[1]);
	memcpy(name, command[1].data(), command[1].size() < MAX_NAME_SIZE ? command[1].size() : MAX_NAME_SIZE);
	if (StringCompare::caseInsensitiveEq(command[2], "NODE")) {
		entity = Entity::NODE;
	}
	if (StringCompare::caseInsensitiveEq(command[2], "ELEMENT") || StringCompare::caseInsensitiveEq(command[2], "ELEM") || StringCompare::caseInsensitiveEq(command[2], "ELEMENTS")) {
		entity = Entity::ELEMENT;
	}
	if (StringCompare::caseInsensitiveEq(Parser::strip(command[2]), "KP")) {
		entity = Entity::KP;
	}
	if (StringCompare::caseInsensitiveEq(Parser::strip(command[2]), "LINE")) {
		entity = Entity::LINE;
	}
	if (StringCompare::caseInsensitiveEq(Parser::strip(command[2]), "AREA")) {
		entity = Entity::AREA;
	}
	if (StringCompare::caseInsensitiveEq(Parser::strip(command[2]), "VOLU")) {
		entity = Entity::VOLU;
	}
	if (entity == Entity::UNKNOWN) {
		error(commandLine);
	}
	NUMITEMS = std::stol(command[3]);

	valueSize = std::stoi(description[0]);
	valueLength = std::stoi(description[1]);

	lineSize = valueSize * valueLength + lineEndSize;

	esint lastLineSize = 0;
	if ((NUMITEMS % valueSize)) {
		lastLineSize = (NUMITEMS % valueSize) * valueLength + lineEndSize;
	}

	WorkbenchParser::fillIndices(begin,
			begin + commandLine.size() + descriptionLine.size(),
			begin + commandLine.size() + descriptionLine.size() + (NUMITEMS / valueSize) * lineSize + lastLineSize);
	return *this;
}

bool CMBlock::readData(std::vector<esint> &indices)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	esint size = (last - first) / lineSize;

	std::vector<esint> tdistribution = tarray<esint>::distribute(threads, size);
	std::vector<std::vector<esint> > tindices(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1];) {
			for (esint n = 0; n < valueSize; ++n) {
				memcpy(value.data(), data, valueLength);
				data += valueLength;
				esint index = atol(value.data());
				if (index != -1 && (t || tindices[t].size() || index > 0)) {
					tindices[t].push_back(index);
				}
			}
			data += lineEndSize;
		}
		if (t == threads - 1) {
			if (lRank == info::mpi::rank) {
				auto data = first + lineSize * tdistribution[t + 1];
				for (esint n = 0; n < NUMITEMS % valueSize; ++n) {
					memcpy(value.data(), data, valueLength);
					data += valueLength;
					esint index = atol(value.data());
					if (index != -1) {
						tindices[t].push_back(index);
					}
				}
			} else {
				// we need to check if the next line starts with a negative value (that indicate array)
				memcpy(value.data(), first + lineSize * tdistribution[t + 1], valueLength);
				esint index = atol(value.data());
				if (index < 0) {
					tindices[t].push_back(index);
				}
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		for (size_t i = 0; i < tindices[t].size(); ++i) {
			if (tindices[t][i] < 0) {
				for (esint ii = indices.back() + 1; ii < -tindices[t][i]; ++ii) {
					indices.push_back(ii);
				}
			} else {
				indices.push_back(tindices[t][i] - 1);
			}
		}
	}

	utils::sortAndRemoveDuplicates(indices);

	return true;
}



