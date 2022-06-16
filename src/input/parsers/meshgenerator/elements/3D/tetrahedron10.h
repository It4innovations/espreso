
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_

#include "quadraticvolume.h"

namespace espreso {

struct Tetrahedron10Generator: public QuadraticVolumeGenerator {

	Tetrahedron10Generator();

	void pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const;
	void pushFace(std::vector<esint> &elements, std::vector<esint> &esize, std::vector<int> &etype, const std::vector<esint> &indices, CubeFace face) const;

	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeEdge edge) const {}
	void pushNodes(std::vector<esint> &nodes, const std::vector<esint> &indices, CubeFace face) const
	{
		pushTriangleNodes(nodes, indices, face);
	}
};

}


#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_3D_TETRAHEDRON10_H_ */
