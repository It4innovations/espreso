
#ifndef SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_
#define SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_

namespace espreso {
namespace input {

enum class CubeFace {
	X_0,
	X_1,
	Y_0,
	Y_1,
	Z_0,
	Z_1,
	NONE
};

enum class CubeEdge {
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


inline std::ostream& operator<<(std::ostream& os, const CubeFace& face)
{
	switch (face) {
		case CubeFace::X_0: return os << "X == 0";
		case CubeFace::X_1: return os << "X == 1";
		case CubeFace::Y_0: return os << "Y == 0";
		case CubeFace::Y_1: return os << "Y == 1";
		case CubeFace::Z_0: return os << "Z == 0";
		case CubeFace::Z_1: return os << "Z == 1";
		default: return os;
	}
}

inline std::ostream& operator<<(std::ostream& os, const CubeEdge& edge)
{
	switch (edge) {
		case CubeEdge::X_0_Y_0: return os << "X == 0 Y == 0";
		case CubeEdge::X_0_Y_1: return os << "X == 0 Y == 1";
		case CubeEdge::X_1_Y_0: return os << "X == 1 Y == 0";
		case CubeEdge::X_1_Y_1: return os << "X == 1 Y == 1";
		case CubeEdge::X_0_Z_0: return os << "X == 0 Z == 0";
		case CubeEdge::X_0_Z_1: return os << "X == 0 Z == 1";
		case CubeEdge::X_1_Z_0: return os << "X == 1 Z == 0";
		case CubeEdge::X_1_Z_1: return os << "X == 1 Z == 1";
		case CubeEdge::Y_0_Z_0: return os << "Y == 0 Z == 0";
		case CubeEdge::Y_0_Z_1: return os << "Y == 0 Z == 1";
		case CubeEdge::Y_1_Z_0: return os << "Y == 1 Z == 0";
		case CubeEdge::Y_1_Z_1: return os << "Y == 1 Z == 1";
		default: return os;
	}
}

}
}



#endif /* SRC_INPUT_MESHGENERATOR_ELEMENTS_ELEMENT_H_ */
