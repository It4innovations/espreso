
#ifndef SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_
#define SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_

#include "keywords.h"

#include <string>
#include <vector>

namespace espreso {

struct InputFilePack;
struct OrderedMeshDatabase;

class EnsightGeometry {
friend class EnsightVariables;

public:
	EnsightGeometry(InputFilePack &geofile, OrderedMeshDatabase &database);

	void scan();
	void parse();

protected:
	InputFilePack &_geofile;
	OrderedMeshDatabase &_database;

	EnsightKeywords _keywords;
};
}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_ */
