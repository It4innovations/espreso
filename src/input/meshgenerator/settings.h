
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "selections/interval.h"
#include "esbasis.h"
#include "esmesh.h"

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

	Settings(const Configuration &configuration, size_t index, size_t size, std::string prefix="");
	Settings(size_t index, size_t size, std::string prefix="");

	std::string prefix;
	std::vector<Parameter> parameters;

	size_t index;
	size_t size;
	size_t clusterOffset;

	std::map<std::string, double> material1;
	std::map<std::string, double> material2;

	std::map<std::string, Interval> nodes;
	std::map<std::string, Interval> faces;
	std::map<std::string, Interval> edges;
	std::map<std::string, Interval> elements;

	std::map<std::string, std::map<std::string, std::string> > properties;

	bool useMetis;
};

}
}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */

