
#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "esmesh.h"
#include "settings.h"
#include "elements/elements.h"

namespace permoncube {

class Generator {

public:
	virtual void mesh(mesh::Mesh &mesh, const size_t cluster[]) = 0;

	virtual ~Generator() { };

protected:
	Generator(permoncube::Settings settings): _settings(settings) { };

	permoncube::Settings _settings;
};

template<class TElement>
class ElementGenerator: public Generator {

public:
	static void globalNodesCount(const Settings &settings, size_t nodes[]);
	static void clusterNodesCount(const Settings &settings, size_t nodes[]);
	static size_t clusterElementsCount(const Settings &settings);

	ElementGenerator(permoncube::Settings &settings): Generator(settings) { };

	void mesh(mesh::Mesh &mesh, const size_t cluster[]);

	~ElementGenerator() { };
};

}

#include "generator.hpp"


#endif /* GENERATOR_H_ */
