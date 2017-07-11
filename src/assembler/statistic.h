
#ifndef SRC_ASSEMBLER_STATISTIC_H_
#define SRC_ASSEMBLER_STATISTIC_H_

#include <cstddef>
#include <vector>

#include "step.h"

namespace espreso {

class Mesh;
class Region;
enum class ElementType;
enum class Property;

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

	Statistic(ElementType eType, const Mesh &mesh, const std::vector<std::vector<double> > &data, const std::vector<Property> &properties);

	void compute(const Step &step);
	double get(const Region* region, size_t offset, StatisticalData statistics);
	double getMagnitude(const Region* region, StatisticalData statistics)
	{
		return get(region, _dataSize, statistics);
	}

private:
	void computeNodes();
	void computeElements();

	// region x offset x data
	std::vector<std::vector<std::vector<double> > > _results;
	Operation _operation;
	ElementType _eType;
	size_t _dataSize;
	bool _computed;
	Step _step;

	const Mesh &_mesh;
	const std::vector<std::vector<double> > &_data;
	std::vector<Region*> _selection;
};

}



#endif /* SRC_ASSEMBLER_STATISTIC_H_ */
