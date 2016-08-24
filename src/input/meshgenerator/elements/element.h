
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_

namespace espreso {
namespace input {

enum class CubeFaces {
	X_0,
	X_1,
	Y_0,
	Y_1,
	Z_0,
	Z_1,
	NONE
};

enum class CubeEdges {
	X_0_Y_0,
	X_0_Y_1,
	X_1_Y_0,
	X_1_Y_1,
	X_0_Z_0,
	X_0_Z_1,
	X_1_Z_0,
	X_1_Z_1,
	Y_0_Z_0,
	Y_0_Z_1,
	Y_1_Z_0,
	Y_1_Z_1,
	NONE
};

}
}



#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_ */
