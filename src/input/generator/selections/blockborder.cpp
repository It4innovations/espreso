

#include "blockborder.h"
#include "../primitives/blocksetting.h"

#include "esbasis.h"

namespace espreso {
namespace input {

BlockBorder::BlockBorder(const std::string &interval)
{
	std::string str = interval;

	std::vector<std::string> values;
	while (str.size()) {
		size_t begin = str.find_first_of("<(");
		size_t end   = str.find_first_of(">)");
		if (end < begin) {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << str << "'. Wrong brackets.";
		}
		values.push_back(str.substr(begin, end + 1));
		if (end + 1 == str.size()) {
			break;
		}
		str.erase(0, end + 1);
	}

	bool _excludeStart[3], _excludeEnd[3];
	double _start[3], _end[3];
	for (size_t i = 0; i < values.size() && i < 3; i++) {
		_excludeStart[i] = _excludeEnd[i] = false;
		std::vector<std::string> bounds = Parser::split(Parser::strip(values[i]), ",;");
		if (bounds.size() != 2) {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Illegal delimiter.";
		}
		if (bounds[0][0] == '(') {
			_excludeStart[i] = true;
		} else if (bounds[0][0] != '<') {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Unknown bracer '" << bounds[0][0] << "'.";
		}
		if (bounds[1][bounds[1].size() - 1] == ')') {
			_excludeEnd[i] = true;
		} else if (bounds[1][bounds[1].size() - 1] != '>') {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Unknown bracer '" << bounds[1][bounds[1].size() - 1] << "'.";
		}

		bounds[0].erase(bounds[0].begin());
		bounds[1].erase(bounds[1].end() - 1);

		std::stringstream ss1(bounds[0]);
		ss1 >> _start[i];

		std::stringstream ss2(bounds[1]);
		ss2 >> _end[i];
	}

	excludeStart = Triple<bool>(_excludeStart[0], _excludeStart[1], _excludeStart[2]);
	excludeEnd = Triple<bool>(_excludeEnd[0], _excludeEnd[1], _excludeEnd[2]);
	start = Triple<double>(_start[0], _start[1], _start[2]);
	end = Triple<double>(_end[0], _end[1], _end[2]);
}

size_t BlockBorder::dimension() const
{
	size_t dimension = 3;
	if (start.x == end.x) { dimension--; }
	if (start.y == end.y) { dimension--; }
	if (start.z == end.z) { dimension--; }
	return dimension;
}

bool BlockBorder::intersect(const BlockSetting &block) const
{
	auto eq = [&] (const double &a, const double &b) {
		return a < b + epsilon && a > b - epsilon ? true : false;
	};
	auto lo = [&] (const double &a, const double &b) {
		return a < b - epsilon ? true : false;
	};

	if (lo(end.x, block.start.x) || (eq(end.x, block.start.x) && excludeEnd.x)) {
		return false;
	}

	if (lo(block.end.x, start.x) || (eq(block.end.x, start.x) && excludeStart.x)) {
		return false;
	}

	if (lo(end.y, block.start.y) || (eq(end.y, block.start.y) && excludeEnd.y)) {
		return false;
	}

	if (lo(block.end.y, start.y) || (eq(block.end.y, start.y) && excludeStart.y)) {
		return false;
	}

	if (lo(end.z, block.start.z) || (eq(end.z, block.start.z) && excludeEnd.z)) {
		return false;
	}

	if (lo(block.end.z, start.z) || (eq(block.end.z, start.z) && excludeStart.z)) {
		return false;
	}

	return true;
}

CubeEdge BlockBorder::getEdge(const BlockSetting &block) const
{
	auto eq = [&] (const double &a, const double &b) {
		return a < b + epsilon && a > b - epsilon ? true : false;
	};
	auto lo = [&] (const double &a, const double &b) {
		return a < b + epsilon ? true : false;
	};

	size_t fixed_x = (start.x == end.x) ? 1 : 0;
	size_t fixed_y = (start.y == end.y) ? 1 : 0;
	size_t fixed_z = (start.z == end.z) ? 1 : 0;

	if (fixed_x + fixed_y + fixed_z < 2) {
		return CubeEdge::NONE;
	}

	if (fixed_z && !eq(block.start.z, start.z) && !eq(block.end.z, end.z)) {
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Warning: interval 'z' coordinate is inside a block. Region is skipped!";
	}
	if (fixed_y && !eq(block.start.y, start.y) && !eq(block.end.y, end.y)) {
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Warning: interval 'y' coordinate is inside a block. Region is skipped!";
	}
	if (fixed_x && !eq(block.start.x, start.x) && !eq(block.end.x, end.x)) {
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Warning: interval 'x' coordinate is inside a block. Region is skipped!";
	}

	if (fixed_z && fixed_y && lo(start.x, block.start.x) && lo(block.start.x, end.x)) {
		if (eq(block.start.z, start.z) && eq(block.start.y, start.y)) {
			return CubeEdge::Y_0_Z_0;
		}
		if (eq(block.start.z, start.z) && eq(block.end.y, end.y)) {
			return CubeEdge::Y_1_Z_0;
		}
		if (eq(block.end.z, end.z) && eq(block.start.y, start.y)) {
			return CubeEdge::Y_0_Z_1;
		}
		if (eq(block.end.z, end.z) && eq(block.end.y, end.y)) {
			return CubeEdge::Y_1_Z_1;
		}
	}

	if (fixed_z && fixed_x && lo(start.y, block.start.y) && lo(block.start.y, end.y)) {
		if (eq(block.start.z, start.z) && eq(block.start.x, start.x)) {
			return CubeEdge::X_0_Z_0;
		}
		if (eq(block.start.z, start.z) && eq(block.end.x, end.x)) {
			return CubeEdge::X_1_Z_0;
		}
		if (eq(block.end.z, end.z) && eq(block.start.x, start.x)) {
			return CubeEdge::X_0_Z_1;
		}
		if (eq(block.end.z, end.z) && eq(block.end.x, end.x)) {
			return CubeEdge::X_1_Z_1;
		}
	}

	if (fixed_y && fixed_x && lo(start.z, block.start.z) && lo(block.start.z, end.z)) {
		if (eq(block.start.y, start.y) && eq(block.start.x, start.x)) {
			return CubeEdge::X_0_Y_0;
		}
		if (eq(block.start.y, start.y) && eq(block.end.x, end.x)) {
			return CubeEdge::X_1_Y_0;
		}
		if (eq(block.end.y, end.y) && eq(block.start.x, start.x)) {
			return CubeEdge::X_0_Y_1;
		}
		if (eq(block.end.y, end.y) && eq(block.end.x, end.x)) {
			return CubeEdge::X_1_Y_1;
		}
	}

	return CubeEdge::NONE;
}

CubeFace BlockBorder::getFace(const BlockSetting &block) const
{
	auto eq = [&] (const double &a, const double &b) {
		return a < b + epsilon && a > b - epsilon ? true : false;
	};
	auto lo = [&] (const double &a, const double &b) {
		return a < b + epsilon ? true : false;
	};

	size_t fixed_x = (start.x == end.x) ? 1 : 0;
	size_t fixed_y = (start.y == end.y) ? 1 : 0;
	size_t fixed_z = (start.z == end.z) ? 1 : 0;

	if (fixed_x + fixed_y + fixed_z != 1) {
		return CubeFace::NONE;
	}

	if (fixed_z && !eq(block.start.z, start.z) && !eq(block.end.z, end.z)) {
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Warning: interval 'z' coordinate is inside a block. Region is skipped!";
	}
	if (fixed_y && !eq(block.start.y, start.y) && !eq(block.end.y, end.y)) {
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Warning: interval 'y' coordinate is inside a block. Region is skipped!";
	}
	if (fixed_x && !eq(block.start.x, start.x) && !eq(block.end.x, end.x)) {
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Warning: interval 'x' coordinate is inside a block. Region is skipped!";
	}

	if (fixed_x && lo(start.y, block.start.y) && lo(block.start.y, end.y) && lo(start.z, block.start.z) && lo(block.start.z, end.z)) {
		if (eq(block.start.x, start.x)) {
			return CubeFace::X_0;
		}
		if (eq(block.end.x, end.x)) {
			return CubeFace::X_1;
		}
	}

	if (fixed_y && lo(start.x, block.start.x) && lo(block.start.x, end.x) && lo(start.z, block.start.z) && lo(block.start.z, end.z)) {
		if (eq(block.start.y, start.y)) {
			return CubeFace::Y_0;
		}
		if (eq(block.end.y, end.y)) {
			return CubeFace::Y_1;
		}
	}

	if (fixed_z && lo(start.x, block.start.x) && lo(block.start.x, end.x) && lo(start.y, block.start.y) && lo(block.start.y, end.y)) {
		if (eq(block.start.z, start.z)) {
			return CubeFace::Z_0;
		}
		if (eq(block.end.z, end.z)) {
			return CubeFace::Z_1;
		}
	}

	return CubeFace::NONE;
}

}
}
