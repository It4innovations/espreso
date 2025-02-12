
#include "operator.h"

#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

void ElementOperatorBuilder::now()
{
    double start;
    if (Operator::print > 1) {
        start = eslog::time();
    }

    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
            for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
                apply(i);
            }
        }
    }

    if (Operator::print > 1) {
        printf("EVALUATE %s: %fs\n", name, eslog::time() - start);
    }
}

void BoundaryOperatorBuilder::now()
{
    double start;
    if (Operator::print > 1) {
        start = eslog::time();
    }

    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
            for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
                if (info::mesh->boundaryRegions[r]->dimension) {
                    for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
                        apply(r, i);
                    }
                }
            }
        }
    }

    if (Operator::print > 1) {
        printf("EVALUATE %s: %fs\n", name, eslog::time() - start);
    }
}
