
#include "vector.sparse.feti.h"
#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"

#include <vector>

using namespace espreso;

VectorSparseFETI::VectorSparseFETI()
{

}

VectorSparseFETI::~VectorSparseFETI()
{

}

VectorSparseFETI* VectorSparseFETI::copy()
{
	return new VectorSparseFETI();
}

Vector* VectorSparseFETI::create()
{
	return new VectorSparse();
}

void VectorSparseFETI::filter(const DataDecomposition *other, esint nindices, esint *indices)
{
//	DataDecomposition::shallowCopy(other);
//	DataDecompositionFiltered::filter(other, nindices, indices);
//
//	std::vector<esint> dsize(domains);
//	for (auto di = dmap->datatarray().cbegin(); di != dmap->datatarray().cend(); ++di) {
//		if (ismy(di->domain)) {
//			++dsize[di->domain - doffset];
//		}
//	}
//	std::vector<esint> ssize(domains);
//	for (auto di = sparsemap->datatarray().cbegin(); di != sparsemap->datatarray().cend(); ++di) {
//		if (ismy(di->domain)) {
//			++ssize[di->domain - doffset];
//		}
//	}
//	for (esint i = 0, d = distribution[rank]; d < distribution[rank + 1]; ++d, ++i) {
//		at(i)->resize(dsize[i], ssize[d]);
//		esint offset = 0, *vals = at(i)->indices;
//		for (auto di = sparsemap->datatarray().begin(); di != sparsemap->datatarray().end(); ++di) {
//			if (di->domain == d) {
//				vals[offset] = di->index;
//				di->index = offset;
//				++offset;
//			}
//		}
//	}
}

VectorsSparseFETI::VectorsSparseFETI()
{
	initVectors(0);
}


VectorsSparseFETI::~VectorsSparseFETI()
{

}

VectorsSparseFETI* VectorsSparseFETI::copy()
{
	return new VectorsSparseFETI();
}

void VectorsSparseFETI::initDomains(DataDecomposition::DUPLICATION duplications, esint domains)
{
	holder()->initDomains(duplications, domains);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->initDomains(duplications, domains);
	}
}

void VectorsSparseFETI::resizeDomain(esint domain, esint size, esint nnz)
{
	holder()->at(domain)->resize(size * nvectors, nnz * nvectors);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->at(domain)->shallowCopyFromHolder(holder()->at(domain), n, nvectors);
	}
}

void VectorsSparseFETI::filter(const DataDecomposition *other, esint nindices, esint *indices)
{
	holder()->filter(other, nindices, indices);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(holder(), n, nvectors);
	}
}

Vector* VectorsSparseFETI::create()
{
	return new VectorSparseFETI();
}

