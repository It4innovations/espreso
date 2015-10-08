
#ifndef INPUT_MESHGENERATOR_CUBE_SETTINGS_H_
#define INPUT_MESHGENERATOR_CUBE_SETTINGS_H_

#include <map>

#include "../configuration/configuration.h"
#include "esmesh.h"

namespace esinput {

struct CubeSettings {

	enum Faces {
		FRONT,
		REAR,
		LEFT,
		RIGHT,
		TOP,
		BOTTOM
	};

	CubeSettings(int argc, char** argv);

	static std::vector<Description> description;

	size_t clusters[3];
	size_t subdomainsInCluster[3];
	size_t elementsInSubdomain[3];
	double problemLength[3];

	std::vector<std::map<size_t, double> > boundaryCondition;
	std::vector<std::map<size_t, bool> > fillCondition;
};

inline std::ostream& operator<<(std::ostream& os, const CubeSettings &s)
{
	os << "clusters: " << s.clusters[0] << " : " << s.clusters[1] << " : " << s.clusters[2] << "\n";
	os << "subdomainsInCluster: " << s.subdomainsInCluster[0] << " : " << s.subdomainsInCluster[1] << " : " << s.subdomainsInCluster[2] << "\n";
	os << "elementsInSubdomain: " << s.elementsInSubdomain[0] << " : " << s.elementsInSubdomain[1] << " : " << s.elementsInSubdomain[2] << "\n";
	os << "cube length: " << s.problemLength[0] << " : " << s.problemLength[1] << " : " << s.problemLength[2] << "\n";
	return os;
}

}

#endif /* INPUT_MESHGENERATOR_CUBE_SETTINGS_H_ */
