
#include "constrains.h"
#include "equalityconstrains.h"
#include "fixedwall.h"

#include "esinfo/ecfinfo.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void Constrains<T>::set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	EqualityConstrains<T>::set(step, feti, dirichlet);
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		FixedWall<T>::set(step, feti); break;
	}
}

template <typename T>
void Constrains<T>::update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	EqualityConstrains<T>::update(step, feti, dirichlet);
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		FixedWall<T>::update(step, feti); break;
	}
}

template struct Constrains<double>;

}

