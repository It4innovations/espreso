
#ifndef SRC_INPUT_PARSERS_ENSIGHT_PARSER_VARIABLES_H_
#define SRC_INPUT_PARSERS_ENSIGHT_PARSER_VARIABLES_H_

#include "casefile.h"
#include "geometry.h"

namespace espreso {

struct InputFilePack;
struct OrderedMeshDatabase;

class EnsightVariables {
public:
	EnsightVariables(const EnsightCasefile &casefile, const EnsightGeometry &geofile, InputFilePack &variables);

	void scan();
	void parse();

protected:
	const EnsightCasefile &_casefile;
	const EnsightGeometry &_geofile;
	InputFilePack &_variables;

	EnsightKeywords _keywords;
};
}



#endif /* SRC_INPUT_PARSERS_ENSIGHT_PARSER_VARIABLES_H_ */
