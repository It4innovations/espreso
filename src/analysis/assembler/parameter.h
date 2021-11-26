
#ifndef SRC_PHYSICS_ASSEMBLER_PARAMETER_H_
#define SRC_PHYSICS_ASSEMBLER_PARAMETER_H_

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "mesh/store/nameddata.h"

#include <vector>

namespace espreso {

class Evaluator;
class ECFExpression;

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

template <class Settings>
struct ParameterSettings {
	std::vector<const Settings*> settings;

	ParameterSettings(int intervals): settings(intervals, NULL) {}
};

struct ParameterData {
	PerElementSize size;
	serializededata<esint, double>* data;

	ParameterData(PerElementSize mask, int intervals);

	void setConstness(bool constness);
	void setUpdate(int value);

	virtual int increment(int interval) const =0;
	virtual int increment(PerElementSize size, int interval) const =0;
	virtual void resize(double init = .0) =0;
	virtual void resizeAligned(size_t alignment, double init = .0) =0;
	virtual ~ParameterData();

	std::string name;
	std::vector<int> isconst, update;
};

struct ExternalElementValue {
	int dimension;
	std::vector<Evaluator*> evaluator;

	ExternalElementValue(ParameterData &parameter);

protected:
	ExternalElementValue(int dimension): dimension(dimension) {}
};

struct ExternalElementGPsValue: public ExternalElementValue { using ExternalElementValue::ExternalElementValue; };
struct ExternalElementNodesValue: public ExternalElementValue { using ExternalElementValue::ExternalElementValue; };

struct ElementParameterData: public ParameterData {
	ElementParameterData(PerElementSize mask);

	static int intervals();
	int increment(int interval) const;
	int increment(PerElementSize size, int interval) const;
	void resize(double init = .0);
	void resizeAligned(size_t alignment, double init = .0);
};

template<int mask>
struct ElementParameter: public ElementParameterData {
	ElementParameter(): ElementParameterData(static_cast<ElementSizeMask>(mask)) { }
};

template<int mask>
struct ElementGPsExternalParameter: public ElementParameter<mask> {
	ExternalElementGPsValue externalValues;

	ElementGPsExternalParameter(): externalValues(*this)
	{
		std::fill(this->update.begin(), this->update.end(), -1);
	}

	bool isSet(size_t interval)
	{
		return this->update[interval] != -1;
	}
};

template<int mask>
struct ElementNodesExternalParameter: public ElementParameter<mask> {
	ExternalElementNodesValue externalValues;

	ElementNodesExternalParameter(): externalValues(*this)
	{
		std::fill(this->update.begin(), this->update.end(), false);
	}
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

struct ExternalBoundaryValue: public ExternalElementValue { // the same as ExternalElementValue but array is array of regions
	ExternalBoundaryValue(BoundaryParameterPack &parameter);
};

template<int mask>
struct BoundaryParameter: public BoundaryParameterPack {
	BoundaryParameter(): BoundaryParameterPack(static_cast<ElementSizeMask>(mask)) { }
};

template<int mask>
struct BoundaryExternalParameter: public BoundaryParameter<mask> {
	ExternalBoundaryValue externalValues;

	BoundaryExternalParameter(): externalValues(*this) { }
};

template <class Settings>
struct BoundaryParameterSettings: public ParameterSettings<Settings> {
	BoundaryParameterSettings(int region): ParameterSettings<Settings>(1), region(region) {}

	int region;
};

template <class Settings>
struct BoundarySettingsPack {
	std::vector<BoundaryParameterSettings<Settings> > regions;

	BoundarySettingsPack()
	{
		regions.reserve(BoundaryParameterData::regions());
		for (int r = 0; r < BoundaryParameterData::regions(); ++r) {
			regions.emplace_back(BoundaryParameterSettings<Settings>(r));
		}
	}
};

template<class Settings>
struct BoundarySettings: public BoundarySettingsPack<Settings> {
	BoundarySettings(): BoundarySettingsPack<Settings>() { }
};

struct InputParameterIterator {
	const int inc;
	const double * __restrict__ data;

	InputParameterIterator(const ParameterData &info, esint interval)
	: inc(info.isconst[interval] ? 0 : info.increment(info.size, interval)), data((info.data->begin() + interval)->data()) {}
	InputParameterIterator(const ParameterData &info, esint interval, PerElementSize size)
	: inc(info.isconst[interval] ? 0 : info.increment(size, interval)), data((info.data->begin() + interval)->data()) {}

	inline InputParameterIterator& operator++() { data += inc; return *this; }
	inline InputParameterIterator& operator+=(const int rhs) { data += rhs*inc; return *this; }
	inline InputParameterIterator& operator+=(const size_t rhs) { data += rhs*inc; return *this; }
	inline const double& operator[](esint i) const { return data[i]; }
};

struct OutputParameterIterator {
	const int inc;
	double * __restrict__ data;

	OutputParameterIterator(ParameterData &info, esint interval)
	: inc(info.isconst[interval] ? 0 : info.increment(info.size, interval)), data((info.data->begin() + interval)->data()) { }
	OutputParameterIterator(ParameterData &info, esint interval, PerElementSize size)
	: inc(info.isconst[interval] ? 0 : info.increment(size, interval)), data((info.data->begin() + interval)->data()) { }

	inline OutputParameterIterator& operator++() { data += inc; return *this; }
	inline OutputParameterIterator& operator+=(const int rhs) { data += rhs*inc; return *this; }
	inline OutputParameterIterator& operator+=(const size_t rhs) { data += rhs*inc; return *this; }
	inline double& operator[](esint i) { return data[i]; }
	inline const double& operator[](esint i) const { return data[i]; }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_PARAMETER_H_ */
