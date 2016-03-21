
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "configuration/configuration.h"

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

struct Settings {

	Settings(const Options &options, size_t index, size_t size);
	Settings(size_t index, size_t size);

	static std::vector<Description> description;

	size_t index;
	size_t size;

	eslocal elementType;
	eslocal shape;

	bool useMetis;
};

inline std::ostream& operator<<(std::ostream& os, const Settings &s)
{
	std::vector<std::string> shapes({ "CUBE", "SPHERE" });
	std::vector<std::string> eTypes({ "HEXA8", "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" });

	os << "clusters: " << s.size << "\n";
	os << "generated shape: " << shapes[s.shape] << "\n";
	os << "type of the element: " << eTypes[s.elementType] << "\n";
	os << "partition: " << (s.useMetis ? "by METIS" : "regular") << "\n";
	return os;
}

}
}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */

