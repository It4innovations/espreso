
#include "faceList.h"

#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"

using namespace espreso;

void OpenFOAMFaceList::parse(ivector<Element::CODE> &type, ivector<esint> &enodes)
{
	type.resize(distribution.back().offset + distribution.back().size);
	std::vector<std::vector<esint> > nodes(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<esint> tnodes;
		tnodes.reserve(4 * distribution[t].size); // more space for general polygons?

		char* next;
		const char *c = distribution[t].c;
		for (size_t i = distribution[t].offset; i < distribution[t].offset + distribution[t].size; ++i) {
			int size = strtol(c, &next, 10); c = next + 1; // skip "("
			switch (size) {
			case 3: type[i] = Element::CODE::TRIANGLE3; break;
			case 4: type[i] = Element::CODE::SQUARE4; break;
			default: type[i] = Element::CODE::POLYGON; tnodes.push_back(size);
			}
			for (esint n = 0; n < size; ++n) {
				tnodes.push_back(strtol(c, &next, 10)); c = next;
			}
			c += 2; // skip ")\n"
		}

		nodes[t].swap(tnodes);
	}

	std::vector<esint> ncount(info::env::OMP_NUM_THREADS + 1);
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		ncount[t] = nodes[t].size();
	}
	utils::sizesToOffsets(ncount);
	enodes.resize(ncount.back());

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::copy(nodes[t].begin(), nodes[t].end(), enodes.begin() + ncount[t]);
	}
}
