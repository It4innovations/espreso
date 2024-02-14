
#include "regularization.h"
#include "regularization.heattransfer.h"
#include "regularization.elasticity.h"
#include "regularization.empty.h"

#include "esinfo/ecfinfo.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void Regularization<T>::set(const step::Step &step, FETI<T> &feti)
{
	feti.R1.resize(feti.K.size());
	feti.R2.resize(feti.K.size());
	feti.RegMat.resize(feti.K.size());

	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
		switch (info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1).type) {
		case HeatTransferLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationHeatTransfer<T>::set(feti); break;
		case HeatTransferLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::set(feti); break;
		} break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
		switch (info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1).type) {
		case StructuralMechanicsLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationElasticity<T>::set(feti); break;
		case StructuralMechanicsLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::set(feti); break;
		case StructuralMechanicsLoadStepConfiguration::TYPE::HARMONIC:     break;
		} break;
	}
}

template <typename T>
void Regularization<T>::update(const step::Step &step, FETI<T> &feti)
{
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
		switch (info::ecf->heat_transfer.load_steps_settings.at(step.loadstep + 1).type) {
		case HeatTransferLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationHeatTransfer<T>::update(feti); break;
		case HeatTransferLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::update(feti); break;
		} break;

	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
		switch (info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1).type) {
		case StructuralMechanicsLoadStepConfiguration::TYPE::STEADY_STATE: RegularizationElasticity<T>::update(feti); break;
		case StructuralMechanicsLoadStepConfiguration::TYPE::TRANSIENT:    RegularizationEmpty<T>::update(feti); break;
		case StructuralMechanicsLoadStepConfiguration::TYPE::HARMONIC:     break;
		} break;
	}
}

template struct Regularization<double>;

}
