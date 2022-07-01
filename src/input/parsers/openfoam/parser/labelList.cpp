
#include "labelList.h"

#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

using namespace espreso;

void OpenFOAMLabelList::parse(ivector<esint> &list)
{
	list.resize(distribution.back().offset + distribution.back().size);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		char *next;
		const char *c = distribution[t].c;
		for (size_t i = distribution[t].offset; i < distribution[t].offset + distribution[t].size; ++i) {
			list[i] = strtol(c, &next, 10); c = next;
			c += 1; // skip "\n"
		}
	}
}

FoamFileHeader OpenFOAMLabelList::load(const std::string &file, ivector<esint> &list, esint offset)
{
	std::ifstream is(file);
	FoamFileHeader header(is);

	esint size;
	is >> size;
	list.reserve(list.size() + size);
	is.ignore(256, '(');
	if (header.format == FoamFileHeader::Format::ASCII) {
		for (esint i = 0; i < size; ++i) {
			esint number;
			is >> number;
			list.push_back(number + offset);
		}
	} else {
		eslog::error("implement OpenFOAM binary reader.\n");
	}
	return header;
}

