
#include "boundary.h"

#include "input/parsers/openfoam/openfoam.h"

#include "basis/utilities/parser.h"

using namespace espreso;

bool OpenFOAMBoundary::readData(OpenFOAMData &data)
{
	esint nFaces = 0, startFace = 0;
	std::string name, parameter;

	const char *c = begin - 3;
	while (*c != '\n') { c--; } // go before number of boundaries

	int n = readInteger(c);
	while (*c++ != '('); // skip '('

	for (int i = 0; i < n; i++) {
		name = readString(c);

		while (*c++ != '{');
		while (*c != '}') {
			parameter = readString(c);
			if (StringCompare::caseSensitiveEq(parameter, "nFaces")) {
				nFaces = readInteger(c);
			}
			if (StringCompare::caseSensitiveEq(parameter, "startFace")) {
				startFace = readInteger(c);
			}
			while (*c == ';' || isEmpty(c)) { ++c; }
		}
		++c;
		while (isEmpty(c)) { ++c; }

		auto &indices = data.eregions[name];
		for (size_t f = 0; f < data.fIDs.size(); f++) {
			if (startFace <= data.fIDs[f] && data.fIDs[f] < startFace + nFaces) {
				indices.push_back(data.fIDs[f]);
			}
		}
	}

	return true;
}


