
#ifndef SRC_OUTPUT_STORE_H_
#define SRC_OUTPUT_STORE_H_

#include "esmesh.h"

namespace espreso {
namespace output {

class Store {

public:
	enum class ElementType {
		NODES,
		EDGES,
		FACES,
		ELEMENTS
	};

	virtual void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shringCluster) = 0;

	virtual void storeGeometry(size_t timeStep = -1) = 0;
	virtual void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType) = 0;
	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType) = 0;
	virtual void finalize() =0;

	virtual ~Store() {};

protected:
	Store(const Mesh &mesh, const std::string &path, double shrinkSubdomain = config::output::SUBDOMAINS_SHRINK_RATIO, double shringCluster = config::output::CLUSTERS_SHRINK_RATIO)
	:_mesh(mesh), _path(path), _shrinkSubdomain(shrinkSubdomain), _shringCluster(shringCluster) {};

	const Mesh &_mesh;
	std::string _path;
	double _shrinkSubdomain;
	double _shringCluster;
};

}
}



#endif /* SRC_OUTPUT_STORE_H_ */
