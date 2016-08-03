
#ifndef INPUT_MESHGENERATOR_SETTINGS_H_
#define INPUT_MESHGENERATOR_SETTINGS_H_

#include "esbasis.h"
#include "esmesh.h"

namespace espreso {
namespace input {

struct Settings {

	Settings(const Configuration &configuration, size_t index, size_t size);
	Settings(size_t index, size_t size);

	std::vector<Parameter> parameters;

	size_t index;
	size_t size;

	std::vector<Material> materials;

	bool useMetis;

protected:
	void defaultSettings();
};

inline std::ostream& operator<<(std::ostream& os, const Settings &s)
{
       os << "clusters: " << s.size << "\n";

       os << "material 1:\n";
       os << "\tdensity: " << s.materials[0].density << "\n";
       os << "\tyoung's modulus: " << s.materials[0].youngModulus << "\n";
       os << "\tpoisson's ratio: " << s.materials[0].poissonRatio << "\n";

       os << "material 2:\n";
       os << "\tdensity: " << s.materials[1].density << "\n";
       os << "\tyoung's modulus: " << s.materials[1].youngModulus << "\n";
       os << "\tpoisson's ratio: " << s.materials[1].poissonRatio << "\n";

       os << "partition: " << (s.useMetis ? "by METIS" : "regular") << "\n";
       return os;
}

}
}



#endif /* INPUT_MESHGENERATOR_SETTINGS_H_ */

