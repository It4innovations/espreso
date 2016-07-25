
#include "interval.h"

#include <vector>
#include <string>

using namespace espreso;

std::ostream& espreso::operator<<(std::ostream& os, const espreso::Interval& interval)
{
	for (size_t i = 0; i < 3; i++) {
		os << (interval.excludeStart[i] ? "(" : "<");
		os << interval.start[i];
		os << ",";
		os << interval.end[i];
		os << (interval.excludeEnd[i] ? ")" : ">");
		os << " ";
	}
    return os;
}

std::istream& espreso::operator>>(std::istream& is, espreso::Interval& interval)
{
	std::string str;
	getline(is, str);

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

	for (size_t i = 0; i < values.size() && i < 3; i++) {
		std::vector<std::string> bounds = Parser::split(Parser::strip(values[i]), ",;");
		if (bounds.size() != 2) {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Illegal delimiter.";
		}
		if (bounds[0][0] == '(') {
			interval.excludeStart[i] = true;
		} else if (bounds[0][0] != '<') {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Unknown bracer '" << bounds[0][0] << "'.";
		}
		if (bounds[1][bounds[1].size() - 1] == ')') {
			interval.excludeEnd[i] = true;
		} else if (bounds[1][bounds[1].size() - 1] != '>') {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Unknown bracer '" << bounds[1][bounds[1].size() - 1] << "'.";
		}

		bounds[0].erase(bounds[0].begin());
		bounds[1].erase(bounds[1].end() - 1);

		std::stringstream ss1(bounds[0]);
		ss1 >> interval.start[i];

		std::stringstream ss2(bounds[1]);
		ss2 >> interval.end[i];
	}

	return is;
}
