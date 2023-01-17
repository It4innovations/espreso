
#ifndef SRC_BASIS_CONTAINERS_EDATA_H_
#define SRC_BASIS_CONTAINERS_EDATA_H_

#include <cstddef>
#include <cstring>

namespace espreso {

template <typename TData>
class edata {

	template <typename TEBoundaries, typename TEData> friend class serializededata;
public:
	size_t size() const          { return this->_end - this->_begin; }

	inline TData& operator[] (size_t n) { return *(this->_begin + n); }
	inline TData& at(size_t n)          { return *(this->_begin + n); }
	inline TData* data()                { return   this->_begin; }
	inline TData& back()                { return *(this->_end - 1); }
	inline TData& front()               { return * this->_begin; }
	inline TData* begin()               { return   this->_begin; }
	inline TData* end()                 { return   this->_end; }

	inline const TData& operator[] (size_t n) const { return *(this->_begin + n); }
	inline const TData* data()                const { return   this->_begin; }

	edata(TData *data, size_t begin, size_t end): _begin(data + begin), _end(data + end) {}
private:
	edata(TData *begin, TData *end): _begin(begin), _end(end) {}
	TData *_begin;
	TData *_end;
};

template <typename TData>
bool operator==(const edata<TData> &d1, const edata<TData> &d2)
{
	if (d1.size() != d2.size()) {
		return false;
	}
	return memcmp(d1.data(), d2.data(), sizeof(TData) * d1.size()) == 0;
}

template <typename TData>
bool operator!=(const edata<TData> &d1, const edata<TData> &d2)
{
	return !(d1 == d2);
}

template <typename TData>
bool operator<(const edata<TData> &d1, const edata<TData> &d2)
{
	for (size_t i = 0; i < d1.size() && i < d2.size(); ++i) {
		if (d1[i] != d2[i]) {
			return d1[i] < d2[i];
		}
	}
	return d1.size() < d2.size();
}

}



#endif /* SRC_BASIS_CONTAINERS_EDATA_H_ */
