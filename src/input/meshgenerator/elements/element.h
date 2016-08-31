
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


inline std::ostream& operator<<(std::ostream& os, const CubeFaces& face)
{
	switch (face) {
		case CubeFaces::X_0: return os << "X == 0";
		case CubeFaces::X_1: return os << "X == 1";
		case CubeFaces::Y_0: return os << "Y == 0";
		case CubeFaces::Y_1: return os << "Y == 1";
		case CubeFaces::Z_0: return os << "Z == 0";
		case CubeFaces::Z_1: return os << "Z == 1";
		default: return os;
	}
}

inline std::ostream& operator<<(std::ostream& os, const CubeEdges& edge)
{
	switch (edge) {
		case CubeEdges::X_0_Y_0: return os << "X == 0 Y == 0";
		case CubeEdges::X_0_Y_1: return os << "X == 0 Y == 1";
		case CubeEdges::X_1_Y_0: return os << "X == 1 Y == 0";
		case CubeEdges::X_1_Y_1: return os << "X == 1 Y == 1";
		case CubeEdges::X_0_Z_0: return os << "X == 0 Z == 0";
		case CubeEdges::X_0_Z_1: return os << "X == 0 Z == 1";
		case CubeEdges::X_1_Z_0: return os << "X == 1 Z == 0";
		case CubeEdges::X_1_Z_1: return os << "X == 1 Z == 1";
		case CubeEdges::Y_0_Z_0: return os << "Y == 0 Z == 0";
		case CubeEdges::Y_0_Z_1: return os << "Y == 0 Z == 1";
		case CubeEdges::Y_1_Z_0: return os << "Y == 1 Z == 0";
		case CubeEdges::Y_1_Z_1: return os << "Y == 1 Z == 1";
		default: return os;
	}
}

}
}



#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_ */
