

#include "sets.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "esinfo/envinfo.h"

#include <fstream>
#include <glob.h>
#include <algorithm>

using namespace espreso;

OpenFOAMSet::OpenFOAMSet()
{
	memset(name, '\0', MAX_NAME_SIZE);
	type = SetType::CELL_SET;
}

OpenFOAMSet::OpenFOAMSet(const std::string &name, SetType type)
: type(type)
{
	memset(this->name, '\0', MAX_NAME_SIZE);
	memcpy(this->name, name.data(), name.size() < MAX_NAME_SIZE ? name.size() : MAX_NAME_SIZE);
}

void OpenFOAMSets::inspect(const std::string &path, std::vector<OpenFOAMSet> &sets)
{
	glob_t glob_result;
	glob(path.c_str(), GLOB_TILDE,NULL, &glob_result);
	for(size_t i = 0; i < glob_result.gl_pathc; ++i){
		std::string parameter, value, name;
		OpenFOAMSet::SetType type = OpenFOAMSet::SetType::CELL_SET;
		std::ifstream is(glob_result.gl_pathv[i]);
		while (is.get() != '{');
		while (is.peek() != '}') {
			is >> parameter >> value;
			while (is.get() != '\n');
			if (StringCompare::caseInsensitiveEq(parameter, "class")) {
				if (StringCompare::caseInsensitiveEq(value, "cellSet;")) {
					type = OpenFOAMSet::SetType::CELL_SET;
				}
				if (StringCompare::caseInsensitiveEq(value, "faceSet;")) {
					type = OpenFOAMSet::SetType::FACE_SET;
				}
				if (StringCompare::caseInsensitiveEq(value, "pointSet;")) {
					type = OpenFOAMSet::SetType::POINT_SET;
				}
			}
			if (StringCompare::caseInsensitiveEq(parameter, "object")) {
				name = std::string(value.begin(), value.end() - 1);
			}
		}
		is.close();
		sets.push_back(OpenFOAMSet(name, type));
	}
}

bool OpenFOAMSets::readData(OpenFOAMSet &set, std::vector<esint> &indices)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);
	std::vector<std::vector<esint> > data(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;

		const char *c = begin + tdistribution[t];
		if (c > begin) {
			while (c < end && *(c - 1) != '\n') { ++c; }
		}
		while (c < begin + tdistribution[t + 1]) {
			tdata.push_back(readInteger(c));
			c += 1; // skip '\n'
		}

		data[t].swap(tdata);
	}

	for (size_t t = 0; t < threads; t++) {
		indices.insert(indices.end(), data[t].begin(), data[t].end());
	}
	std::sort(indices.begin(), indices.end());

	return true;
}











