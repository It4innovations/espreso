
#include "constrains.h"
#include "equalityconstrains.h"
#include "fixedwall.h"

#include "esinfo/ecfinfo.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void Constrains<T>::set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	eq.set(step, feti, dirichlet);
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER: break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
		fw.set(step, feti, dirichlet); break;
	}
}

template <typename T>
void Constrains<T>::update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
	eq.update(step, feti, dirichlet);
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER: break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
		fw.update(step, feti, dirichlet); break;
	}
}

template struct Constrains<double>;

}

