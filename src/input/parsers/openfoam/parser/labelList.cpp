
#include "labelList.h"

#include "esinfo/envinfo.h"

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
