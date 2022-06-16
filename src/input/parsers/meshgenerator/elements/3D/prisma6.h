
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_

#include "linearvolume.h"

namespace espreso {

struct Prisma6Generator: public LinearVolumeGenerator {

	Prisma6Generator();

	void pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const;
	void pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const;

	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const {}
	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const;
};

}

#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_PRISMA6_H_ */
