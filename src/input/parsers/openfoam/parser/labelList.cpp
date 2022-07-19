
#include "labelList.h"

#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

using namespace espreso;

template <int T> static esint _parse(std::ifstream &is);

template <> esint _parse<4>(std::ifstream &is)
{
	int number;
	is.read(reinterpret_cast<char*>(&number), 4);
	return number;
}

template <> esint _parse<8>(std::ifstream &is)
{
	long number;
	is.read(reinterpret_cast<char*>(&number), 8);
	return number;
}

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
	FoamFileHeader header; header.read(is);

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
		for (esint i = 0; i < size; ++i) {
			switch (header.label) {
			case 4: list.push_back(_parse<4>(is) + offset); break;
			case 8: list.push_back(_parse<8>(is) + offset); break;
			}
		}
	}
	return header;
}

