
#ifndef SRC_INPUT_FORMATS_ENSIGHT_PARSER_CASEFILE_H_
#define SRC_INPUT_FORMATS_ENSIGHT_PARSER_CASEFILE_H_

#include <string>
#include <vector>

namespace espreso {

class EnsightCasefile {
public:
	enum class Type {
		UNKNOWN,
		ENSIGHT_GOLD
	};

	struct Variable {
		enum class Type {
			NODE,
			ELEMENT
		};

		int dimension, timeset;
		Type type;
		std::string name, path;
		int skip;

		Variable(int dimension, Type type, int timeset, const std::string &name, const std::string &path)
		: dimension(dimension), timeset(timeset), type(type), name(name), path(path), skip(false)
		{

		}
	};

	struct TimeSet {
		int tstart, tend, tinc;
		std::vector<double> values;
	};

	EnsightCasefile(const std::string &path);

	void parse();

	std::string path;
	Type type;
	std::string geometry;
	std::vector<Variable> variables;
	std::vector<TimeSet> timesets;
};

}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_PARSER_CASEFILE_H_ */
