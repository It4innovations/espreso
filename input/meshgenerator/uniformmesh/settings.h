
#ifndef INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_
#define INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_

#include "../configuration/configuration.h"
#include "esmesh.h"

namespace esinput {

struct UniformSettings {

	UniformSettings(int argc, char** argv);

	static std::vector<Description> description;

	size_t subdomainsInCluster[3];
	size_t elementsInSubdomain[3];
	size_t cornerCount;
	bool corners;
	bool edges;
	bool faces;
};

inline std::ostream& operator<<(std::ostream& os, const UniformSettings &s)
{
	os << "subdomainsInCluster: " << s.subdomainsInCluster[0] << " : " << s.subdomainsInCluster[1] << " : " << s.subdomainsInCluster[2] << "\n";
	os << "elementsInSubdomain: " << s.elementsInSubdomain[0] << " : " << s.elementsInSubdomain[1] << " : " << s.elementsInSubdomain[2] << "\n";
	return os;
}

}



#endif /* INPUT_MESHGENERATOR_UNIFORMMESH_SETTINGS_H_ */
