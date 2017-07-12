
#ifndef SRC_CONFIGURATION_MATERIAL_HOLDER_H_
#define SRC_CONFIGURATION_MATERIAL_HOLDER_H_

#include "../configuration.hpp"
#include "parameters.h"

namespace espreso {

template <MATERIAL_PARAMETER TParameter>
struct MaterialParam {

	MaterialParam() {}
	MaterialParam(const std::string &&value) : value(value) {}
	std::string value;
};

template <MATERIAL_PARAMETER TParameter>
struct ValueHolder<MaterialParam<TParameter> >: public Parameter {
	MaterialParam<TParameter> &value;

	ValueHolder(const std::string &name, const std::string &description, MaterialParam<TParameter> &value, MaterialParam<TParameter> defaultValue, const std::string &type)
	: Parameter(name, description, "*"), value(value) {};

	bool set(const std::string &value)
	{
		this->value.value = value;
		return true;
	}

	std::string get() const
	{
		return value.value;
	}

	size_t index() const { return (size_t)TParameter; }

	std::string XMLAttributeType() const { return "string"; }
};

}



#endif /* SRC_CONFIGURATION_MATERIAL_HOLDER_H_ */
