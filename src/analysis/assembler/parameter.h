
#ifndef SRC_PHYSICS_ASSEMBLER_PARAMETER_H_
#define SRC_PHYSICS_ASSEMBLER_PARAMETER_H_

#include "basis/containers/serializededata.h"

#include <vector>
#include <algorithm>

namespace espreso {

enum ElementSizeMask: int {
	ndim   = 1 << 6,
	edim   = 1 << 12,
	enodes = 1 << 18,
	egps   = 1 << 24,
};

inline constexpr ElementSizeMask operator*(const int &i1, const ElementSizeMask &i2) { return (ElementSizeMask)(i1 + i2); }
inline constexpr ElementSizeMask operator*(const ElementSizeMask &i1, const int &i2) { return (ElementSizeMask)(i1 + i2); }
inline constexpr ElementSizeMask operator*(const ElementSizeMask &i1, const ElementSizeMask &i2) { return (ElementSizeMask)(i1 + i2); }

struct PerElementSize {
	int n, ndimension, edimension, node, gp;

	PerElementSize(int mask): n(63 & mask ? 63 & mask : 1), ndimension(63 & (mask >> 6)),  edimension(63 & (mask >> 12)), node(63 & (mask >> 18)), gp(63 & (mask >> 24)) {}
	PerElementSize(ElementSizeMask mask): n(63 & mask ? 63 & mask : 1), ndimension(63 & (mask >> 6)), edimension(63 & (mask >> 12)), node(63 & (mask >> 18)), gp(63 & (mask >> 24)) {}

	bool operator==(const PerElementSize &other) const { return n == other.n && ndimension == other.ndimension && edimension == other.edimension && node == other.node && gp == other.gp; }
	bool operator!=(const PerElementSize &other) const { return !(*this == other); }
};

struct ParameterData {
	PerElementSize size;
	serializededata<esint, double>* data;

	ParameterData(PerElementSize mask, int intervals);

	void setConstness(bool constness);

	virtual int increment(int interval) const =0;
	virtual int increment(PerElementSize size, int interval) const =0;
	virtual void resize(double init = .0) =0;
	virtual ~ParameterData();

	std::vector<int> isconst;
};

struct ElementParameterData: public ParameterData {
	ElementParameterData(PerElementSize mask);

	static int intervals();
	int increment(int interval) const;
	int increment(PerElementSize size, int interval) const;
	void resize(double init = .0);
};

template<int mask>
struct ElementParameter: public ElementParameterData {
	ElementParameter(): ElementParameterData(static_cast<ElementSizeMask>(mask)) { }
};

struct BoundaryParameterData: public ParameterData {
	BoundaryParameterData(int region, PerElementSize mask);

	static int regions();
	static int intervals(int region);
	int increment(int interval) const;
	int increment(PerElementSize size, int interval) const;
	void resize(double init = .0);
	void resizeAligned(size_t alignment, double init = .0);

	bool isSet()
	{
		return data != nullptr;
	}

	int region;
};

struct BoundaryParameterPack {
	std::vector<BoundaryParameterData> regions;

	BoundaryParameterPack(PerElementSize mask);

	bool isSet(size_t region)
	{
		return regions[region].isSet();
	}

	PerElementSize size;
};

template<int mask>
struct BoundaryParameter: public BoundaryParameterPack {
	BoundaryParameter(): BoundaryParameterPack(static_cast<ElementSizeMask>(mask)) { }
};


struct InputParameterIterator {
	const int inc;
	const double * __restrict__ data;

	InputParameterIterator(const ParameterData &info, esint interval)
	: inc(info.isconst[interval] ? 0 : info.increment(info.size, interval)), data(info.data ? (info.data->begin() + interval)->data() : nullptr) {}
	InputParameterIterator(const ParameterData &info, esint interval, PerElementSize size)
	: inc(info.isconst[interval] ? 0 : info.increment(size, interval)), data(info.data ? (info.data->begin() + interval)->data() : nullptr) {}

	inline InputParameterIterator& operator++() { data += inc; return *this; }
	inline InputParameterIterator& operator+=(const int rhs) { data += rhs*inc; return *this; }
	inline InputParameterIterator& operator+=(const size_t rhs) { data += rhs*inc; return *this; }
	inline const double& __restrict__ operator[](esint i) const { return data[i]; }
};

struct OutputParameterIterator {
	const int inc;
	double * __restrict__ data;

	OutputParameterIterator(ParameterData &info, esint interval)
	: inc(info.isconst[interval] ? 0 : info.increment(info.size, interval)), data((info.data->begin() + interval)->data()) { }
	OutputParameterIterator(ParameterData &info, esint interval, PerElementSize size)
	: inc(info.isconst[interval] ? 0 : info.increment(size, interval)), data((info.data->begin() + interval)->data()) { }

	inline OutputParameterIterator& operator++() { data += inc; return *this; }
	inline OutputParameterIterator& operator+=(const int rhs) { data += rhs * inc; return *this; }
	inline OutputParameterIterator& operator+=(const size_t rhs) { data += rhs * inc; return *this; }
	inline double& __restrict__ operator[](esint i) { return data[i]; }
	inline const double& __restrict__ operator[](esint i) const { return data[i]; }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_PARAMETER_H_ */
