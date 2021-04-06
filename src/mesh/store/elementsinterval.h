
#ifndef SRC_MESH_STORE_ELEMENTSINTERVAL_H_
#define SRC_MESH_STORE_ELEMENTSINTERVAL_H_

namespace espreso {

struct ElementsInterval {
	esint begin, end;
	esint domain;
	int code;

	ElementsInterval(): begin(0), end(0), domain(-1), code(-1) {}
	ElementsInterval(esint begin, esint end): begin(begin), end(end), domain(-1), code(-1) {}
	ElementsInterval(esint begin, esint end, esint domain, int code) : begin(begin), end(end), domain(domain), code(code) {}

	bool operator==(const ElementsInterval &other) const { return begin == other.begin && end == other.end && domain == other.domain && code == other.code; }
};

}


#endif /* SRC_MESH_STORE_ELEMENTSINTERVAL_H_ */
