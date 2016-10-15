
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "selections/interval.h"
#include "esbasis.h"
#include "esmesh.h"

namespace espreso {
namespace input {

struct Settings {

	Settings(const ArgsConfiguration &configuration, size_t index, size_t size, std::string prefix="");
	Settings(size_t index, size_t size, std::string prefix="");

	std::string prefix;
	std::vector<Parameter> parameters;

	size_t index;
	size_t size;
	size_t clusterOffset;

	std::map<std::string, std::string> material1;
	std::map<std::string, std::string> material2;

	std::map<std::string, Interval> nodes;
	std::map<std::string, Interval> faces;
	std::map<std::string, Interval> edges;
	std::map<std::string, Interval> elements;

	std::map<std::string, std::map<std::string, std::string> > properties;

	bool useMetis;

protected:
	void defaultSettings();
};

}
}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */

