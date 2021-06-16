
#ifndef SRC_MATH2_UTILS_UTILS_FETI_H_
#define SRC_MATH2_UTILS_UTILS_FETI_H_

#include <vector>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;

struct DI { esint domain, index; };
inline bool operator==(const DI &left, const DI &right) { return left.domain == right.domain && left.index == right.index; }
inline bool operator!=(const DI &left, const DI &right) { return !(left == right); }
inline bool operator<(const DI &left, const DI &right)  { return left.domain == right.domain ? left.index < right.index : left.domain < right.domain; }

class DomainDecomposition
{
public:
	enum class DUPLICATION {
		SPLIT,
		DUPLICATE,
		SPLIT_DOMAINS
	};

	bool ismy(int rank, esint domain) const
	{
		return distribution[rank] <= domain && domain < distribution[rank + 1];
	}

	DUPLICATION duplications;
	std::vector<esint> distribution;
	std::vector<int> neighbors;
	serializededata<esint, DI> *map;
};


}

#endif /* SRC_MATH2_UTILS_UTILS_FETI_H_ */
