

#include "blockborder.h"
#include "../primitives/blocksetting.h"
#include "../generator.h"

#include "../../../basis/logging/logging.h"
#include "../../../basis/utilities/parser.h"

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
	start = (Triple<double>(_start[0] / Generator::precision, _start[1] / Generator::precision, _start[2] / Generator::precision)).round();
	end = (Triple<double>(_end[0] / Generator::precision, _end[1] / Generator::precision, _end[2] / Generator::precision)).round();
}

size_t BlockBorder::dimension() const
{
	if (start.x == end.x && (excludeStart.x || excludeEnd.x)) { return -1; }
	if (start.y == end.y && (excludeStart.y || excludeEnd.y)) { return -1; }
	if (start.z == end.z && (excludeStart.z || excludeEnd.z)) { return -1; }

	size_t dimension = 3;
	if (start.x == end.x) { dimension--; }
	if (start.y == end.y) { dimension--; }
	if (start.z == end.z) { dimension--; }
	return dimension;
}

BlockBorder BlockBorder::intersect(const BlockSetting &block) const
{
	BlockBorder intersection(*this);

	if (intersection.start.x < block.start.x) { intersection.start.x = block.start.x; intersection.excludeStart.x = false; }
	if (intersection.start.y < block.start.y) { intersection.start.y = block.start.y; intersection.excludeStart.y = false; }
	if (intersection.start.z < block.start.z) { intersection.start.z = block.start.z; intersection.excludeStart.z = false; }

	if (intersection.end.x > block.end.x) { intersection.end.x = block.end.x; intersection.excludeEnd.x = false; }
	if (intersection.end.y > block.end.y) { intersection.end.y = block.end.y; intersection.excludeEnd.y = false; }
	if (intersection.end.z > block.end.z) { intersection.end.z = block.end.z; intersection.excludeEnd.z = false; }

	if (intersection.end.x < intersection.start.x) { intersection.end.x = intersection.start.x; intersection.excludeStart.x = true; }
	if (intersection.end.y < intersection.start.y) { intersection.end.y = intersection.start.y; intersection.excludeStart.y = true; }
	if (intersection.end.z < intersection.start.z) { intersection.end.z = intersection.start.z; intersection.excludeStart.z = true; }

	return intersection;
}

CubeEdge BlockBorder::getEdge(const BlockSetting &block) const
{
	auto eq = [&] (const esglobal &a, const esglobal &b) {
		return a == b;
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

	if (fixed_z && fixed_y) {
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

	if (fixed_z && fixed_x) {
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

	if (fixed_y && fixed_x) {
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
		return a == b;
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

	if (fixed_x) {
		if (eq(block.start.x, start.x)) {
			return CubeFace::X_0;
		}
		if (eq(block.end.x, end.x)) {
			return CubeFace::X_1;
		}
	}

	if (fixed_y) {
		if (eq(block.start.y, start.y)) {
			return CubeFace::Y_0;
		}
		if (eq(block.end.y, end.y)) {
			return CubeFace::Y_1;
		}
	}

	if (fixed_z) {
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
