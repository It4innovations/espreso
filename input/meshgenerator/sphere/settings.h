
#ifndef INPUT_MESHGENERATOR_SPHERE_SETTINGS_H_
#define INPUT_MESHGENERATOR_SPHERE_SETTINGS_H_

namespace esinput {

struct SphereSettings {

	SphereSettings() {
		clusters[0] = clusters[1] = clusters[2] = 1;
		subdomainsInCluster[0] = subdomainsInCluster[1] = subdomainsInCluster[2] = 2;
		elementsInSubdomain[0] = elementsInSubdomain[1] = elementsInSubdomain[2] = 5;
		layers = 1;
	};

	size_t clusters[3];
	size_t subdomainsInCluster[3];
	size_t elementsInSubdomain[3];
	size_t layers;
};

inline std::ostream& operator<<(std::ostream& os, const SphereSettings &s)
{
	os << "clusters: " << s.clusters[0] << " : " << s.clusters[1] << " : " << s.clusters[2] << "\n";
	os << "subdomainsInCluster: " << s.subdomainsInCluster[0] << " : " << s.subdomainsInCluster[1] << " : " << s.subdomainsInCluster[2] << "\n";
	os << "elementsInSubdomain: " << s.elementsInSubdomain[0] << " : " << s.elementsInSubdomain[1] << " : " << s.elementsInSubdomain[2] << "\n";
	os << "layers: " << s.layers << "\n";
	return os;
}

}


#endif /* INPUT_MESHGENERATOR_SPHERE_SETTINGS_H_ */
