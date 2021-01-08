
#include "operator.h"

#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

void ElementOperatorBuilder::now()
{
	if (Operator::print) printf("BUILD ELEMENTS: %s\n", name());
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				apply(i);
			}
		}
	}
}

void BoundaryOperatorBuilder::now()
{
	if (Operator::print) printf("BUILD BOUNDARY: %s\n", name());
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
						apply(r, i);
					}
				}
			}
		}
	}
}
