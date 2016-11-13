
#include "interval.h"

#include <vector>
#include <string>

namespace espreso {

std::ostream& operator<<(std::ostream& os, const Interval& interval)
{
	for (size_t i = 0; i < 3; i++) {
		os << (interval.excludeStart[i] ? "(" : "<");
		os << interval.start[i] << "," << interval.end[i];
		os << (interval.excludeEnd[i] ? ")" : ">");
	}
	return os;
}

std::istream& operator>>(std::istream& is, Interval& interval)
{
	std::string str;
	getline(is, str);

	if (StringCompare::caseInsensitiveEq(Parser::strip(str), "ALL")) {
		return is;
	}
	interval._all = false;

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

	bool excludeStart, excludeEnd;
	for (size_t i = 0; i < values.size() && i < 3; i++) {
		excludeStart = excludeEnd = false;
		std::vector<std::string> bounds = Parser::split(Parser::strip(values[i]), ",;");
		if (bounds.size() != 2) {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Illegal delimiter.";
		}
		if (bounds[0][0] == '(') {
			excludeStart = true;
		} else if (bounds[0][0] != '<') {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Unknown bracer '" << bounds[0][0] << "'.";
		}
		if (bounds[1][bounds[1].size() - 1] == ')') {
			excludeEnd = true;
		} else if (bounds[1][bounds[1].size() - 1] != '>') {
			ESINFO(GLOBAL_ERROR) << "Cannot parse interval '" << values[i] << "'. Unknown bracer '" << bounds[1][bounds[1].size() - 1] << "'.";
		}

		bounds[0].erase(bounds[0].begin());
		bounds[1].erase(bounds[1].end() - 1);

		std::stringstream ss1(bounds[0]);
		ss1 >> interval.start[i];
		interval.excludeStart[i] = excludeStart;

		std::stringstream ss2(bounds[1]);
		ss2 >> interval.end[i];
		interval.excludeEnd[i] = excludeEnd;
	}

	return is;
}

}
