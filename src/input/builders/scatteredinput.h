
#ifndef SRC_INPUT_SCATTEREDINPUT_H_
#define SRC_INPUT_SCATTEREDINPUT_H_

#include "input.h"
#include "basis/sfc/hilbertcurve.h"

namespace espreso {

class ScatteredInput: public Input {
public:
	ScatteredInput(MeshBuilder &dMesh);

protected:
	void assignNBuckets();
	void assignEBuckets();

	void clusterize();
	void computeSFCNeighbors();
	void mergeDuplicatedNodes();
	void linkup();
	void exchangeBoundary();

	HilbertCurve<esfloat> _sfc;

	std::vector<esint> _nIDs;
	std::vector<esint> _nBuckets, _eBuckets;

	// distribution across processes
	std::vector<esint> _bucketsBorders;
	std::vector<int> _sfcNeighbors;
};

}


#endif /* SRC_INPUT_SCATTEREDINPUT_H_ */
