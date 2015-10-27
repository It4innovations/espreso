
#ifndef PHYSICS_LINEAR_TEMPERATURE_TEMPERATURE_H_
#define PHYSICS_LINEAR_TEMPERATURE_TEMPERATURE_H_

#include "../linear.h"

namespace physics {

template <MatrixComposer TMatrixComposer>
class Temperature: public Linear<TMatrixComposer> {

public:
	Temperature(const mesh::Mesh &mesh): Linear<TMatrixComposer>(mesh) {};

	void init();
	void pre_solve_update();
	void post_solve_update();
	void solve();
	void finalize();

protected:
	size_t DOFs() { return 1; }
	void inertia(std::vector<double> &inertia)
	{
		inertia.resize(1, 0);
	}
	void C(DenseMatrix &C);
	double CP() { return 1; }
	double rho() { return 7.85e-9; }
};

}

#include "linearelasticity.hpp"




#endif /* PHYSICS_LINEAR_TEMPERATURE_TEMPERATURE_H_ */
