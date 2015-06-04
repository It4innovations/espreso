
#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "esmesh.h"
#include "settings.h"
#include "elements/elements.h"
#include "utils.h"

namespace permoncube {

class Generator {

public:
	virtual void mesh(mesh::Mesh &mesh, const size_t cluster[]) = 0;
	virtual void fixZeroPlanes(
			std::map<eslocal, double>  &dirichlet_x,
			std::map<eslocal, double>  &dirichlet_y,
			std::map<eslocal, double>  &dirichlet_z,
			const size_t cluster[]) = 0;
	virtual void fixBottom(
			std::map<eslocal, double>  &dirichlet_x,
			std::map<eslocal, double>  &dirichlet_y,
			std::map<eslocal, double>  &dirichlet_z,
			const size_t cluster[]) = 0;

	virtual ~Generator() { };

protected:
	Generator(permoncube::Settings settings): _settings(settings) { };

	permoncube::Settings _settings;
};

template<class TElement>
class ElementGenerator: public Generator {

public:
	ElementGenerator(permoncube::Settings &settings): Generator(settings) { };

	void mesh(mesh::Mesh &mesh, const size_t cluster[]);
	void fixZeroPlanes(
			std::map<eslocal, double>  &dirichlet_x,
			std::map<eslocal, double>  &dirichlet_y,
			std::map<eslocal, double>  &dirichlet_z,
			const size_t cluster[]);
	void fixBottom(
			std::map<eslocal, double>  &dirichlet_x,
			std::map<eslocal, double>  &dirichlet_y,
			std::map<eslocal, double>  &dirichlet_z,
			const size_t cluster[]);

	~ElementGenerator()
	{
		TElement::clear();
	}
};

}

#include "generator.hpp"


#endif /* GENERATOR_H_ */
