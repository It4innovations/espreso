
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "esbasis.h"
#include "esmesh.h"

#include "parsers/interval.h"

namespace espreso {
namespace input {

enum GeneratorShape {
	CUBE,
	SPHERE
};

enum ElementType {
	HEXA8,
	HEXA20,
	TETRA4,
	TETRA10,
	PRISMA6,
	PRISMA15,
	PYRAMID5,
	PYRAMID13
};

enum Assembler {
	LinearElasticity,
	Temperature,
	TransientElasticity
};

struct Settings {

	Settings(const Options &options, size_t index, size_t size);
	Settings(size_t index, size_t size);

	std::vector<Description> description;

	size_t index;
	size_t size;

	std::map<std::string, double> material1;
	std::map<std::string, double> material2;

	bool useMetis;
};

}
}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */

