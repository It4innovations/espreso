
#ifndef SRC_CONFIG_MATERIALHOLDER_H_
#define SRC_CONFIG_MATERIALHOLDER_H_

#include "materialparameters.h"
#include "configuration.h"

namespace espreso {

template <MATERIAL_PARAMETER TParameter>
struct MaterialParam {

	MaterialParam() {}
	MaterialParam(const std::string &&value) : value(value) {}
	std::string value;
};

template <MATERIAL_PARAMETER TParameter>
struct ValueHolder<MaterialParam<TParameter> >: public ParameterBase {
	MaterialParam<TParameter> &value;

	ValueHolder(const std::string &name, const std::string &description, MaterialParam<TParameter> &value, MaterialParam<TParameter> defaultValue, const std::string &type)
	: ParameterBase(name, description, "*"), value(value) {};

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
};

}



#endif /* SRC_CONFIG_MATERIALHOLDER_H_ */
