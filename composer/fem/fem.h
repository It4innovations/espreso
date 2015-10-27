
#ifndef COMPOSER_FEM_FEM_H_
#define COMPOSER_FEM_FEM_H_

#include "../composer.h"
#include "esmesh.h"

namespace composer {

template <class TPhysics>
class FEM: public Composer {

public:
	FEM(const mesh::Mesh &mesh): _mesh(mesh) {};

	void element(SparseMatrix &K, size_t subdomain);
	void subdomainBoundary();
	void clusterBoundary();

protected:
	const mesh::Mesh &_mesh;
};

}

#include "fem.hpp"


#endif /* COMPOSER_FEM_FEM_H_ */
