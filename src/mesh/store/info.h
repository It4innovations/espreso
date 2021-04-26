
#ifndef SRC_MESH_STORE_INFO_H_
#define SRC_MESH_STORE_INFO_H_

namespace espreso {

struct UniqueDataInfo {
	esint offset, last, size, totalSize;

	UniqueDataInfo(): offset(0), last(0), size(0), totalSize(0) {}

	inline bool isLocal(const esint &index) const noexcept { return offset <= index && index < last; }
};

struct DuplicateDataInfo: public UniqueDataInfo {
	esint nhalo;

	DuplicateDataInfo(): nhalo(0) {}
};

}

#endif /* SRC_MESH_STORE_INFO_H_ */
