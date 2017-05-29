
#ifndef SRC_ASSEMBLER_STATISTIC_H_
#define SRC_ASSEMBLER_STATISTIC_H_

#include <cstddef>
#include <vector>

namespace espreso {

class Mesh;
class Region;
enum class ElementType;

enum StatisticalData: int {
	MIN     = 1 << 0,
	MAX     = 1 << 1,
	AVERAGE = 1 << 2,
	NORM    = 1 << 3,
	SQUARES = 1 << 4
};

class Statistic {

public:


	// how to treat data distributed to more domains
	enum class Operation: int {
		SUM,
		AVERAGE
	};

	Statistic(StatisticalData statistics, Operation operation, ElementType eType, const Mesh &mesh, const std::vector<std::vector<double> > &data, const std::vector<size_t> &offsets, const std::vector<Region*> &selection);
	Statistic(ElementType eType, const Mesh &mesh, const std::vector<std::vector<double> > &data, size_t dataSize);

	void compute();
	double get(const Region* region, size_t offset, StatisticalData statistics);

private:
	// region x offset x data
	std::vector<std::vector<std::vector<double> > > _results;
	StatisticalData _statistics;
	Operation _operation;
	ElementType _eType;
	size_t _dataSize;
	bool _computed;

	const Mesh &_mesh;
	const std::vector<std::vector<double> > &_data;
	std::vector<size_t> _offsets;
	std::vector<Region*> _selection;
};

inline constexpr StatisticalData  operator& (StatisticalData  d1, StatisticalData  d2) { return StatisticalData( static_cast<int>(d1) & static_cast<int>(d2)); }
inline constexpr StatisticalData  operator| (StatisticalData  d1, StatisticalData  d2) { return StatisticalData( static_cast<int>(d1) | static_cast<int>(d2)); }
inline constexpr StatisticalData  operator^ (StatisticalData  d1, StatisticalData  d2) { return StatisticalData( static_cast<int>(d1) ^ static_cast<int>(d2)); }
inline constexpr StatisticalData  operator~ (StatisticalData  d1)                      { return StatisticalData(~static_cast<int>(d1)                       ); }

inline const     StatisticalData& operator&=(StatisticalData &d1, StatisticalData &d2) { return d1 = d1 & d2; }
inline const     StatisticalData& operator|=(StatisticalData &d1, StatisticalData &d2) { return d1 = d1 | d2; }
inline const     StatisticalData& operator^=(StatisticalData &d1, StatisticalData &d2) { return d1 = d1 ^ d2; }

}



#endif /* SRC_ASSEMBLER_STATISTIC_H_ */
