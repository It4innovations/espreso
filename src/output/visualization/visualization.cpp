
#include "visualization.h"
#include "esinfo/stepinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "mesh/store/nameddata.h"

using namespace espreso;

Visualization::Visualization(const Mesh &mesh)
: ResultStoreBase(mesh)
{
	if (info::ecf->output.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
//		createOutputDirectory();
	}
}

Visualization::~Visualization()
{

}

bool Visualization::isRoot()
{
	return info::mpi::rank == 0;
//	if (info::mpi::rank == 0) {
//		if (step::type == step::TYPE::FTT) {
//			return true;
//		}
//		if (step::duplicate::instances == 1 && info::mpi::grank == 0) {
//			return true;
//		}
//	}
//	return false;
}

bool Visualization::storeStep()
{
	if (step::type == step::TYPE::FTT) {
		return true;
	} else {
		switch (info::ecf->output.results_store_frequency) {
		case OutputConfiguration::STORE_FREQUENCY::NEVER:
			return false;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_SUBSTEP:
			return true;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_SUBSTEP:
			return step::substep % info::ecf->output.results_nth_stepping == 0;
		case OutputConfiguration::STORE_FREQUENCY::LAST_SUBSTEP:
			return step::isLast();

		default:
			return false;
		}
	}
}

bool Visualization::storeData(const NamedData *data)
{
	if (static_cast<int>(data->restriction & step::type) == 0) {
		return false;
	}
	if (data->name.size() == 0) {
		return false;
	}
	return true;
}

Point Visualization::shrink(const Point &p, const Point &ccenter, const Point &dcenter, double cratio, double dratio) {
	Point point = ccenter + (p - ccenter) * cratio;
	point = dcenter + (point - dcenter) * dratio;
	return point;
}
