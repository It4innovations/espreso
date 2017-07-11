
#ifndef SRC_ASSEMBLER_SOLUTION_H_
#define SRC_ASSEMBLER_SOLUTION_H_

#include <cstddef>
#include <vector>
#include <string>

#include "statistic.h"
#include "step.h"

namespace espreso {

enum class Property;
enum class ElementType;


struct Solution {

	Solution(const Mesh &mesh, const std::string &name, ElementType eType, const std::vector<Property> &properties, const std::vector<std::vector<double> > &data);
	Solution(const Mesh &mesh, const std::string &name, ElementType eType, const std::vector<Property> &properties);

	inline double get(Property property, eslocal domain, eslocal index) const
	{
		return data[domain][index * properties.size() + _offset[static_cast<int>(property)]];
	}

	inline double get(size_t propertyOffset, eslocal domain, eslocal index) const
	{
		return data[domain][index * properties.size() + propertyOffset];
	}

	std::vector<std::vector<double> >& innerData() { return _data; }

	bool hasProperty(Property property) const;
	void computeStatisticalData(const Step &step);
	double getStatisticalData(const std::vector<Property> &property, StatisticalData data, const Region *region) const;

	std::string name;
	ElementType eType;
	std::vector<Property> properties;
	const std::vector<std::vector<double> > &data;

protected:
	std::vector<int> _offset;

	// when no data are provided, store it here
	std::vector<std::vector<double> > _data;

	mutable Statistic _statistic;
};

}



#endif /* SRC_ASSEMBLER_SOLUTION_H_ */
