
#include "dataholder.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/fetidatastore.h"

#include <algorithm>

using namespace espreso;

DataHolder::DataHolder()
: decomposition(NULL)
{

}

void DataHolder::assembleB0fromKernels()
{
	std::vector<esint> rindex;

	auto dual = info::mesh->FETIData->domainDual->begin();
	std::vector<esint> rows(info::mesh->elements->nclusters);
	for (esint d1 = 0; d1 < info::mesh->elements->domains.size; ++d1, ++dual) {
		esint cluster = info::mesh->elements->clusters[d1];
		for (auto dit = dual->begin(); dit != dual->end(); ++dit) {
			if (d1 < *dit) {
				rindex.push_back(rows[cluster]);
				if (N1[d1].cols < N1[*dit].cols) {
					rows[cluster] += N1[*dit].cols;
				} else {
					rows[cluster] += N1[d1].cols ? N1[d1].cols : 1;
				}
			} else {
				auto dualbegin = info::mesh->FETIData->domainDual->begin();
				auto it = std::lower_bound((dualbegin + *dit)->begin(), (dualbegin + *dit)->end(), d1);
				rindex.push_back(rindex[it - dualbegin->begin()]);
			}
		}
	}

	dual = info::mesh->FETIData->domainDual->begin();
	auto dmap = decomposition->dmap->begin();
	for (esint n = 0, prev = 0; n < decomposition->nshared; prev = decomposition->shared[n++]) {
		dmap += decomposition->shared[n] - prev;
		for (auto di1 = dmap->begin(); di1 != dmap->end(); ++di1) {
			for (auto di2 = di1 + 1; di2 != dmap->end(); ++di2) {
				if (decomposition->ismy(di1->domain) && decomposition->ismy(di2->domain)) {
					auto it = std::lower_bound(
							(dual + (di1->domain - decomposition->doffset))->begin(),
							(dual + (di1->domain - decomposition->doffset))->end(),
							di2->domain - decomposition->doffset);

					if (it != (dual + (di1->domain - decomposition->doffset))->end() && *it == di2->domain - decomposition->doffset) {
						esint d1, d2, d1index, d2index;
						if (di1->domain < di2->domain) {
							d1 = di1->domain - decomposition->doffset;
							d2 = di2->domain - decomposition->doffset;
							d1index = di1->index;
							d2index = di2->index;
						} else {
							d1 = di2->domain - decomposition->doffset;
							d2 = di1->domain - decomposition->doffset;
							d1index = di2->index;
							d2index = di1->index;
						}
						if (N1[d1].cols) {
							for (esint c = 0; c < N1[d1].cols; ++c) {
								B0[d1].I_row_indices.push_back(rindex[it - dual->begin()] + c + 1);
								B0[d1].J_col_indices.push_back(d1index + 1);
								B0[d1].V_values.push_back(N1[d1].dense_values[N1[d1].rows * c + d1index]);
								B0[d2].I_row_indices.push_back(rindex[it - dual->begin()] + c + 1);
								B0[d2].J_col_indices.push_back(d2index + 1);
								B0[d2].V_values.push_back(-N1[d1].dense_values[N1[d1].rows * c + d1index]);
							}
						} else {
							B0[d1].I_row_indices.push_back(rindex[it - dual->begin()] + 1);
							B0[d1].J_col_indices.push_back(d1index + 1);
							B0[d1].V_values.push_back( (double)(d1 + 1) / info::mesh->elements->domains.size);
							B0[d2].I_row_indices.push_back(rindex[it - dual->begin()] + 1);
							B0[d2].J_col_indices.push_back(d2index + 1);
							B0[d2].V_values.push_back(-(double)(d1 + 1) / info::mesh->elements->domains.size);
						}
					}
				}
			}
		}
	}

	for (esint d = 0; d < info::mesh->elements->domains.size; ++d) {
		B0[d].rows = rows[info::mesh->elements->clusters[d]];
		B0[d].cols = K[d].cols;
		B0[d].nnz = B0[d].V_values.size();
	}
}



