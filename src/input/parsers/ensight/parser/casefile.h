
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

		int dimension, time;
		Type type;
		std::string name, path;

		Variable(int dimension, Type type, int time, const std::string &name, const std::string &path)
		: dimension(dimension), time(time), type(type), name(name), path(path)
		{

		}
	};

	EnsightCasefile(const std::string &path);

	void parse();

	std::string path;
	Type type;
	std::string geometry;
	std::vector<Variable> variables;
	std::vector<std::vector<double> > times;
};

}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_PARSER_CASEFILE_H_ */
