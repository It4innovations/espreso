
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

	virtual void fixZeroPlanes(mesh::Mesh &mesh, const size_t cluster[]) = 0;
	virtual void fixBottom(mesh::Mesh &mesh, const size_t cluster[]) = 0;

	virtual void fillGlobalBoundaries(mesh::Boundaries &boundaries, const size_t cluster[]) = 0;

	virtual void setFixPoints(mesh::Mesh &mesh, const size_t cluster[]) = 0;
	virtual void setCorners(
			mesh::Boundaries &boundaries,
			const size_t cluster[],
			const size_t number[],
			const bool corners,
			const bool edges,
			const bool surface) = 0;

	virtual ~Generator() { };

protected:
	Generator(permoncube::Settings settings): _settings(settings) { };

	permoncube::Settings _settings;
};

template<class TElement>
class ElementGenerator: public Generator {

public:
	ElementGenerator(permoncube::Settings &settings): Generator(settings), e(settings) { };

	void mesh(mesh::Mesh &mesh, const size_t cluster[]);

	void fixZeroPlanes(mesh::Mesh &mesh, const size_t cluster[]);
	void fixBottom(mesh::Mesh &mesh, const size_t cluster[]);

	void fillGlobalBoundaries(mesh::Boundaries &boundaries, const size_t cluster[]);

	void setFixPoints(mesh::Mesh &mesh, const size_t cluster[]);
	void setCorners(
			mesh::Boundaries &boundaries,
			const size_t cluster[],
			const size_t number[],
			const bool corners,
			const bool edges,
			const bool surface);

private:
	TElement e;
};

}

#include "generator.hpp"


#endif /* GENERATOR_H_ */
