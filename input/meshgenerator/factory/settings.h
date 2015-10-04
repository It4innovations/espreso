
#ifndef INPUT_MESHGENERATOR_FACTORY_SETTINGS_H_
#define INPUT_MESHGENERATOR_FACTORY_SETTINGS_H_

#include "../configuration/configuration.h"

namespace esinput {

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

struct FactorySettings {

	FactorySettings(int argc, char** argv);

	static std::vector<Description> description;

	eslocal elementType;
	eslocal shape;
};

inline std::ostream& operator<<(std::ostream& os, const FactorySettings &s)
{
	os << "generated shape: " << s.shape << "\n";
	os << "type of the element: " << s.elementType << "\n";
	return os;
}

}




#endif /* INPUT_MESHGENERATOR_FACTORY_SETTINGS_H_ */
