
#ifndef SRC_CONFIG_HOLDERS_VALUEHOLDER_H_
#define SRC_CONFIG_HOLDERS_VALUEHOLDER_H_

#include "expression.h"
#include "range.h"
#include "config/configuration.h"
#include "basis/utilities/parser.h"
#include <sstream>

namespace espreso {

template <typename Ttype>
struct ECFValueHolder: public ECFValue {
	Ttype &value;

	ECFValueHolder(Ttype &value): value(value) {}

	std::string getValue() const
	{
		std::stringstream ss;
		ss << value;
		return ss.str();
	}

	bool _setValue(const std::string &value)
	{
		std::stringstream ss(value);
		ss >> this->value;
		return ss.eof() && !ss.fail();
	}

	virtual const void* data() const { return &value; }
};

template <>
inline ECFValueHolder<ECFExpression>::ECFValueHolder(ECFExpression &value)
: value(value)
{
	if (value.value.size()) {
		value.createEvaluator();
	}
}

template <>
inline std::string ECFValueHolder<std::string>::getValue() const
{
	return value;
}

template <>
inline std::string ECFValueHolder<ECFExpression>::getValue() const
{
	return value.value;
}

template <>
inline std::string ECFValueHolder<ECFRange>::getValue() const
{
	return value.min + " : " + value.max + " : " + value.step;
}

template <>
inline bool ECFValueHolder<std::string>::_setValue(const std::string &value)
{
	this->value = value;
	return true;
}

template <>
inline bool ECFValueHolder<ECFExpression>::_setValue(const std::string &value)
{
	this->value.value = Parser::uppercase(value);
	this->value.createEvaluator();
	return this->value.evaluator != NULL;
}

template <>
inline bool ECFValueHolder<ECFRange>::_setValue(const std::string &value)
{
	auto range = Parser::split(value, ":");
	if (range.size() != 3) {
		return false;
	}
	this->value.min = range[0];
	this->value.max = range[1];
	this->value.step = range[2];
	return true;
}

template <>
inline bool ECFValueHolder<bool>::_setValue(const std::string &value)
{
	if (value.size() == 0) {
		this->value = true;
		return true;
	} else {
		if (StringCompare::caseInsensitiveEq(value, "FALSE")) {
			this->value = false;
			return true;
		}
		if (StringCompare::caseInsensitiveEq(value, "TRUE")) {
			this->value = true;
			return true;
		}
		std::stringstream ss(value);
		ss >> this->value;
		return ss.eof();
	}
}

template <>
inline std::string ECFValueHolder<bool>::getValue() const
{
	if (value) {
		return "TRUE";
	} else {
		return "FALSE";
	}
}

template <typename Ttype>
struct ECFEnumHolder: public ECFValue {
	Ttype &value;

	ECFEnumHolder(Ttype &value): value(value) {}

	std::string getValue() const {
		if (metadata.datatype.front() == ECFDataType::OPTION) {
			return metadata.options[static_cast<int>(value)].name;
		}
		if (metadata.datatype.front() == ECFDataType::ENUM_FLAGS) {
			std::string flags;
			for (size_t i = 0; i < metadata.options.size(); i++) {
				if (static_cast<int>(value) & static_cast<int>(1 << i)) {
					if (flags.size()) {
						flags += " | ";
					}
					flags += metadata.options[i].name;
				}
			}
			return flags;
		}
		return "";
	}

	bool _setValue(const std::string &value)
	{
		if (metadata.datatype.front() == ECFDataType::ENUM_FLAGS) {
			std::vector<std::string> values = Parser::split(value, "|");
			if (values.size()) {
				this->value = static_cast<Ttype>(0);
				for (size_t i = 0; i < values.size(); i++) {
					values[i] = Parser::strip(values[i]);
					Ttype prev = this->value;
					for (size_t j = 0; j < metadata.options.size(); j++) {
						if (StringCompare::caseInsensitiveEq(values[i], metadata.options[j].name)) {
							this->value = static_cast<Ttype>(static_cast<int>(this->value) | (int)(1 << j));
							break;
						}
					}
					if (this->value == prev) {
						return false;
					}
				}
				return true;
			}
		}
		size_t index = -1;
		for (size_t i = 0; i < metadata.options.size(); i++) {
			if (StringCompare::caseInsensitiveEq(value, metadata.options[i].name)) {
				index = i;
				break;
			}
		}

		if (index == (size_t)-1) {
			if (!ECFValueHolder<size_t>(index)._setValue(value)) {
				return false;
			}
		}

		if (metadata.datatype.front() == ECFDataType::OPTION) {
			this->value = static_cast<Ttype>(index);
			return true;
		}
		if (metadata.datatype.front() == ECFDataType::ENUM_FLAGS) {
			this->value = static_cast<Ttype>(1 << index);
			return true;
		}
		return false;
	}

	virtual const void* data() const { return &value; }

};

}



#endif /* SRC_CONFIG_HOLDERS_VALUEHOLDER_H_ */
