
#ifndef SRC_MESH_STORE_INFO_H_
#define SRC_MESH_STORE_INFO_H_

namespace espreso {

struct DistributionInfo {
	esint offset, size, totalSize;

	DistributionInfo(): offset(0), size(0), totalSize(0) {}
};

struct DistributedDataInfo: public DistributionInfo {
	esint next;

	DistributedDataInfo(): next(0) {}

	inline bool isLocal(const esint &index) const noexcept { return offset <= index && index < next; }
};

struct DuplicatedDataInfo: public DistributedDataInfo {
	esint nhalo;

	DuplicatedDataInfo(): nhalo(0) {}
};

}

#endif /* SRC_MESH_STORE_INFO_H_ */
