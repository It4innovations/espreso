
#ifndef SRC_PHYSICS_ASSEMBLER_PARAMETER_H_
#define SRC_PHYSICS_ASSEMBLER_PARAMETER_H_

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"

namespace espreso {

class ECFExpression;
class NodeData;
struct ExpressionsToElements;
struct ExpressionsToBoundary;

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

	template <class Parameter, class ... Other>
	typename std::enable_if<sizeof...(Other)>::type
	addInputs(const Parameter &p, const Other& ... other)
	{
		addInput(p);
		addInputs(other...);
	}

	template <class Parameter, class ... Other>
	typename std::enable_if<!sizeof...(Other)>::type
	addInputs(const Parameter &p, const Other& ... other)
	{
		addInput(p);
		resize();
	}

	void addInput(const ParameterData &p);
	void addInput(const serializededata<esint, esint>* p);
	void addInput(const serializededata<esint, Point>* p);
	void addInput(const ECFExpression &exp, int interval);
	void addInput(const NodeData* p);

	void addInput(int interval, const serializededata<esint, Point>* p);

	void setConstness(bool constness);

	virtual int increment(int interval) const =0;
	virtual int increment(PerElementSize size, int interval) const =0;
	virtual void resize() =0;
	virtual ~ParameterData() {}

	bool isset;
	std::vector<bool> isconst;
	std::vector<int> version;
};

struct ElementParameterData: public ParameterData {
	ElementParameterData(PerElementSize mask);

	static int intervals();
	int increment(int interval) const;
	int increment(PerElementSize size, int interval) const;
	void resize();
};

template<int mask>
struct ElementParameter: public ElementParameterData {
	ElementParameter(): ElementParameterData(static_cast<ElementSizeMask>(mask)) { }
};

template<int mask>
struct ElementExternalParameter: public ElementParameter<mask> {
	ExpressionsToElements *builder;

	ElementExternalParameter(): builder(NULL) { }
};

struct BoundaryParameterData: public ParameterData {
	BoundaryParameterData(int region, PerElementSize mask);

	static int regions();
	static int intervals(int region);
	int increment(int interval) const;
	int increment(PerElementSize size, int interval) const;
	void resize();

	int region;
};

struct BoundaryParameterPack {
	std::vector<BoundaryParameterData> regions;

	BoundaryParameterPack(PerElementSize mask);

	template <class Parameter>
	void addInput(const Parameter &p)
	{
		for (size_t r = 0; r < regions.size(); ++r) {
			regions[r].addInput(p);
		}
	}

	template <class ... Parameters>
	void addInputs(const Parameters& ...p)
	{
		for (size_t r = 0; r < regions.size(); ++r) {
			regions[r].addInputs(p...);
		}
	}

private:
	PerElementSize _mask;
};

template<int mask>
struct BoundaryParameter: public BoundaryParameterPack {
	BoundaryParameter(): BoundaryParameterPack(static_cast<ElementSizeMask>(mask)) { }
};

template<int mask>
struct BoundaryExternalParameter: public BoundaryParameter<mask> {
	ExpressionsToBoundary *builder;

	BoundaryExternalParameter(): builder(NULL) { }
};

template <class Settings>
struct BoundaryParameterSettings: public ParameterSettings<Settings> {
	BoundaryParameterSettings(int region): ParameterSettings<Settings>(1), region(region), isset(false) {}

	int region, isset;
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
	const double * __restrict data;

	InputParameterIterator(const double * data, int increment): inc(increment), data(data) {}
	InputParameterIterator(const ParameterData &info, esint interval, PerElementSize size)
	: inc(info.isconst[interval] ? 0 : info.increment(size, interval)), data((info.data->begin() + interval)->data()) {}

	inline InputParameterIterator& operator++() { data += inc; return *this; }
	inline const double& operator[](esint i) const { return data[i]; }
};

struct OutputParameterIterator {
	const int inc;
	double * __restrict data;

	OutputParameterIterator(double * data, int increment): inc(increment), data(data) {}
	OutputParameterIterator(ParameterData &info, esint interval, PerElementSize size)
	: inc(info.isconst[interval] ? 0 : info.increment(size, interval)), data((info.data->begin() + interval)->data()) { }

	inline OutputParameterIterator& operator++() { data += inc; return *this; }
	inline double& operator[](esint i) { return data[i]; }
	inline const double& operator[](esint i) const { return data[i]; }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_PARAMETER_H_ */
