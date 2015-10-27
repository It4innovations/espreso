
#ifndef PHYSICS_PHYSICS_H_
#define PHYSICS_PHYSICS_H_

#include "essolver.h"
#include "esmesh.h"
#include "esbem.h"

namespace physics {

enum MatrixComposer {
	FEM,
	BEM,
	ELMER
};

template <MatrixComposer TMatrixComposer>
class Physics {

public:
	virtual void init() = 0;
	virtual void pre_solve_update() = 0;
	virtual void post_solve_update() = 0;
	virtual void solve() = 0;
	virtual void finalize() = 0;

	virtual ~Physics() {};

protected:
	Physics(const mesh::Mesh &mesh);

	const mesh::Mesh &_mesh;
	const mesh::SurfaceMesh _surface;

	bool _verbose;
	TimeEval _timeStatistics;

};

template<>
Physics<FEM>::Physics(const mesh::Mesh &mesh): _mesh(mesh), _surface(mesh.rank(), mesh.size()), _verbose(true) { }

template<>
Physics<BEM>::Physics(const mesh::Mesh &mesh): _mesh(mesh), _surface(mesh), _verbose(true) { }

template<>
Physics<ELMER>::Physics(const mesh::Mesh &mesh): _mesh(mesh), _surface(mesh.rank(), mesh.size()), _verbose(true) { }

}



#endif /* PHYSICS_PHYSICS_H_ */
