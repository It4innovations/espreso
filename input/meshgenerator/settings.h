
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "configuration/configuration.h"

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
	os << "index: " << s.index << "\n";
	os << "size: " << s.size << "\n";
	os << "generated shape: " << s.shape << "\n";
	os << "type of the element: " << s.elementType << "\n";
	os << "use METIS: " << s.useMetis << "\n";
	return os;
}

}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */
