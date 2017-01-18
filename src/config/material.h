
#ifndef SRC_CONFIG_MATERIAL_H_
#define SRC_CONFIG_MATERIAL_H_

#include "materialparameters.h"

namespace espreso {

template <MATERIAL_PARAMETER TParameter>
struct MaterialParam {

	MaterialParam() {}
	MaterialParam(const std::string &&value) : value(value) {}
	std::string value;
};

template <MATERIAL_PARAMETER TParameter>
std::ostream& operator<<(std::ostream& os, const MaterialParam<TParameter>& parameter)
{
	os << parameter.value;
	return os;
}

template <MATERIAL_PARAMETER TParameter>
std::istream& operator>>(std::istream& is, MaterialParam<TParameter>& parameter)
{
	is >> parameter.value;
	return is;
}

struct MaterialParameters { };

}




#endif /* SRC_CONFIG_MATERIAL_H_ */
