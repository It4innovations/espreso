
#ifndef SRC_ASSEMBLER_SOLUTION_H_
#define SRC_ASSEMBLER_SOLUTION_H_

#include <cstddef>
#include <vector>
#include <string>

namespace espreso {

enum class Property;

struct Solution {

	Solution(const std::string &name, const std::vector<Property> &properties, const std::vector<std::vector<double> > &data);
	Solution(const std::string &name, const std::vector<Property> &properties);

	inline double get(Property property, eslocal domain, eslocal index) const
	{
		return data[domain][index * properties + _offset[static_cast<int>(property)]];
	}

	std::string name;
	size_t properties;
	const std::vector<std::vector<double> > &data;
protected:
	std::vector<int> _offset;

	// when no data are provided, store it here
	std::vector<std::vector<double> > _data;
};

}



#endif /* SRC_ASSEMBLER_SOLUTION_H_ */
