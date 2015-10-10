
#ifndef INPUT_MESHGENERATOR_SPHERE_SETTINGS_H_
#define INPUT_MESHGENERATOR_SPHERE_SETTINGS_H_

#include "../configuration/configuration.h"

namespace esinput {

struct SphereSettings {

	SphereSettings(int argc, char** argv);

	static std::vector<Description> description;

	size_t subdomainsInCluster[3];
	size_t elementsInSubdomain[3];
	size_t layers;
	double innerRadius;
	double outerRadius;
	size_t cornerCount;
	bool corners;
	bool edges;
	bool faces;
};

inline std::ostream& operator<<(std::ostream& os, const SphereSettings &s)
{
	os << "subdomainsInCluster: " << s.subdomainsInCluster[0] << " : " << s.subdomainsInCluster[1] << " : " << s.subdomainsInCluster[2] << "\n";
	os << "elementsInSubdomain: " << s.elementsInSubdomain[0] << " : " << s.elementsInSubdomain[1] << " : " << s.elementsInSubdomain[2] << "\n";
	os << "layers: " << s.layers << "\n";
	return os;
}

}


#endif /* INPUT_MESHGENERATOR_SPHERE_SETTINGS_H_ */
