
#ifndef SRC_WRAPPERS_MATH_DOMAININDICES_H_
#define SRC_WRAPPERS_MATH_DOMAININDICES_H_

namespace espreso {

struct DI {
	esint domain, index;
};

inline bool operator==(const DI &left, const DI &right)
{
	return left.domain == right.domain && left.index == right.index;
}

inline bool operator!=(const DI &left, const DI &right)
{
	return !(left == right);
}

inline bool operator<(const DI &left, const DI &right)
{
	return left.domain == right.domain ? left.index < right.index : left.domain < right.domain;
}

}

#endif /* SRC_WRAPPERS_MATH_DOMAININDICES_H_ */
