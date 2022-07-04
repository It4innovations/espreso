
#ifndef SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_
#define SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_

#include "keywords.h"
#include "input/input.h"

#include <string>
#include <vector>

namespace espreso {

struct InputFilePack;
struct OrderedMeshDatabase;

class EnsightGeometry {
friend class EnsightVariables;

public:
	EnsightGeometry(InputFilePack &geofile);

	void scan();
	void parse(OrderedNodes &nodes, OrderedElements &elements, OrderedRegions &regions);

protected:
	InputFilePack &_geofile;

	EnsightKeywords _keywords;
};
}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_PARSER_GEOMETRY_H_ */
