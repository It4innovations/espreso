
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "configuration/configuration.h"

namespace esinput {

struct Settings {

	Settings(int argc, char** argv, size_t index, size_t size);
	Settings(size_t index, size_t size);

	static std::vector<Description> description;

	size_t index;
	size_t size;

	bool useMetis;
};

inline std::ostream& operator<<(std::ostream& os, const Settings &s)
{
	os << "index: " << s.index << "\n";
	os << "size: " << s.size << "\n";
	os << "use METIS: " << s.useMetis << "\n";
	return os;
}

}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */
