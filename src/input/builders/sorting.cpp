
#include "builder.utils.h"

#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"

#include <algorithm>
#include <numeric>

namespace espreso {
namespace builder {

void groupElementTypes(ClusteredMesh &clustered)
{
	std::vector<esint, initless_allocator<esint> > permutation(clustered.eoffsets.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (const esint &i, const esint &j) {
		if (Mesh::edata[(int)clustered.etype[i]].type != Mesh::edata[(int)clustered.etype[j]].type) {
			return (int)Mesh::edata[(int)clustered.etype[i]].type > (int)Mesh::edata[(int)clustered.etype[j]].type;
		} else {
			return clustered.eoffsets[i] < clustered.eoffsets[j];
		}
	});

	utils::permute(clustered.eoffsets, permutation);
	utils::permute(clustered.etype, permutation);
	std::vector<esint> npermutation(clustered.enodes.size());
	for (size_t i = 0, index = 0; i < permutation.size(); i++) {
		for (esint n = 0; n < Mesh::edata[(int)clustered.etype[i]].nodes; ++n, ++index) {
			npermutation[index] = clustered.edist[permutation[i]] + n;
		}
	}
	utils::permute(clustered.enodes, npermutation);

	for (int type = static_cast<int>(Element::TYPE::VOLUME); type > static_cast<int>(Element::TYPE::POINT); --type) {
		clustered.typeDistribution.push_back(std::lower_bound(clustered.etype.begin(), clustered.etype.end(), type, [&] (Element::CODE code, int type) {
			return static_cast<int>(Mesh::edata[(int)code].type) >= type; }) - clustered.etype.begin()
		);
	}
}


}
}
