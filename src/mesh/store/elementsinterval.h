
#ifndef SRC_MESH_STORE_ELEMENTSINTERVAL_H_
#define SRC_MESH_STORE_ELEMENTSINTERVAL_H_

#include <vector>

namespace espreso {

struct ElementsInterval {
	esint begin, end;
	esint domain;
	int code, material, region, einterval;
	std::vector<int> regions; // in the (hopefully) rare case of intersected regions

	ElementsInterval(): begin(0), end(0), domain(-1), code(-1), material(-1), region(-1), einterval(-1) {}
	ElementsInterval(esint begin, esint end): begin(begin), end(end), domain(-1), code(-1), material(-1), region(-1), einterval(-1) {}
//	ElementsInterval(esint begin, esint end, esint domain, int code, int region) : begin(begin), end(end), domain(domain), code(code), region(region) {}

	bool operator==(const ElementsInterval &other) const { return begin == other.begin && end == other.end && domain == other.domain && code == other.code && region == other.region; }
};

}


#endif /* SRC_MESH_STORE_ELEMENTSINTERVAL_H_ */
