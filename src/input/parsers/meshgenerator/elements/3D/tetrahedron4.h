
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_

#include "linearvolume.h"

namespace espreso {

struct Tetrahedron4Generator: public LinearVolumeGenerator {

	Tetrahedron4Generator();

	void pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const;
	void pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const;

	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const {}
	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
	{
		pushTriangleNodes(nodes, indices, face);
	}
};

}


#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON4_H_ */
