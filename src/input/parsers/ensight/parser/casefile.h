
#ifndef SRC_INPUT_FORMATS_ENSIGHT_PARSER_CASEFILE_H_
#define SRC_INPUT_FORMATS_ENSIGHT_PARSER_CASEFILE_H_

#include <string>
#include <vector>

namespace espreso {

class EnsightCasefile {
public:
	enum class Type {
		Ensight_Gold
	};

	EnsightCasefile(const std::string &path);

	void parse();

	std::string path;
	Type type;
	std::string geometry;
	std::vector<std::string> variables;
	std::vector<double> times;
};

}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_PARSER_CASEFILE_H_ */
