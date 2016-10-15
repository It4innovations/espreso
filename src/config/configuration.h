
#ifndef SRC_CONFIG_CONFIGURATION_H_
#define SRC_CONFIG_CONFIGURATION_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "esbasis.h"
#include "reader.h"

#define MAX_VECTOR_CONFIGURATION_SIZE 10

#define OPTIONS(...) __VA_ARGS__
#define OPTION(type, name, description, value, options) type name = DataHolder::create<type>(#name, description, name, value, options, parameters)
#define PARAMETER(type, name, description, value)       type name = DataHolder::create<type>(#name, description, name, value, parameters)
#define SUBCONFIG(type, name)                           type name = Configuration::create<type>(#name, subconfigurations);

namespace espreso {

struct ParameterBase {
	std::string name;
	std::string description;

	ParameterBase(std::string name, std::string description)
	: name(name), description(description) {}

	virtual bool set(const std::string &value) =0;
	virtual std::string get() const =0;
	virtual ~ParameterBase() {};
};

template <typename Tvalue>
struct Option {
	std::string name;
	std::string description;
	Tvalue value;

	Option(std::string name, Tvalue value, std::string description): name(name), description(description), value(value) {}
};



template <typename Ttype>
struct ValueHolder: public ParameterBase {
	Ttype &value;

	ValueHolder(std::string name, std::string description, Ttype &value, Ttype defaultValue)
	: ParameterBase(name, description), value(value) { value = defaultValue; }

	bool set(const std::string &value)
	{
		std::stringstream ss(value);
		ss >> this->value;
		return ss.eof() && !ss.fail();
	}

	std::string get() const
	{
		std::stringstream ss;
		ss << value;
		return ss.str();
	}
};

template <>
inline bool ValueHolder<bool>::set(const std::string &value)
{
	if (value.size() == 0) {
		this->value = true;
		return true;
	} else {
		std::stringstream ss(value);
		ss >> this->value;
		return ss.eof();
	}
}

template <>
inline bool ValueHolder<std::string>::set(const std::string &value)
{
	this->value = value;
	return true;
}

template <>
inline std::string ValueHolder<std::string>::get() const
{
	return value;
}

template <typename Ttype>
struct OptionsHolder: public ParameterBase {
	Ttype &value;
	std::vector<Option<Ttype> > options;

	OptionsHolder(std::string name, std::string description, Ttype &value, Ttype defaultValue, std::vector<Option<Ttype> > options)
	: ParameterBase(name, description), value(value), options(options) { value = defaultValue; }

	bool set(const std::string &value)
	{
		for (size_t i = 0; i < options.size(); i++) {
			if (StringCompare::caseInsensitiveEq(value, options[i].name)) {
				this->value = options[i].value;
				return true;
			}
		}
		std::stringstream ss(value);
		size_t number;
		ss >> number;
		if (!ss.fail() && ss.eof() && number < options.size()) {
			this->value = options[number].value;
			return true;
		}
		return false;
	}

	std::string get() const
	{
		for (size_t i = 0; i < options.size(); i++) {
			if (options[i].value == value) {
				return options[i].name;
			}
		}
		// this code is not reachable
		return "Unrecognized value";
	}
};

struct DataHolder {
	template <typename Ttype>
	static Ttype create(std::string name, std::string description, Ttype &value, Ttype defaultValue, std::map<std::string, ParameterBase*, StringCompare> &parameters)
	{
		parameters[name] = new ValueHolder<Ttype>(name, description, value, defaultValue);
		return value;
	}

	template <typename Ttype>
	static Ttype create(std::string name, std::string description, Ttype &value, Ttype defaultValue, std::vector<Option<Ttype> > options, std::map<std::string, ParameterBase*, StringCompare> &parameters)
	{
		parameters[name] = new OptionsHolder<Ttype>(name, description, value, defaultValue, options);
		return value;
	}
};

struct Configuration {
	std::map<std::string, ParameterBase*, StringCompare> parameters;
	std::map<std::string, Configuration*, StringCompare> subconfigurations;

	static void read(const std::string &file) { Reader::read(file); }
	static void read(int *argc, char ***argv) { Reader::read(argc, argv); }
	static void print() { Reader::print(); }
	static void store() { Reader::store(); }

	template <typename Ttype>
	static Ttype create(const std::string &name, std::map<std::string, Configuration*, StringCompare> &subconfigurations)
	{
		Ttype configuration;
		subconfigurations[name] = &configuration;
		return configuration;
	}

	Configuration& operator[](const std::string &subconfiguration)
	{
		if (subconfigurations.find(subconfiguration) != subconfigurations.end()) {
			return *subconfigurations.find(subconfiguration)->second;
		} else {
			ESINFO(GLOBAL_ERROR) << "Unknown sub-configuration '" << subconfiguration << "'";
			exit(EXIT_FAILURE);
		}
	}

	bool set(const std::string &parameter, const std::string &value)
	{
		if (parameters.find(parameter) != parameters.end()) {
			return parameters.find(parameter)->second->set(value);
		} else {
			ESINFO(GLOBAL_ERROR) << "Unknown parameter '" << parameter << "'";
			return false;
		}
	}

	~Configuration()
	{
		for (auto it = parameters.begin(); it != parameters.end(); ++it) {
			delete it->second;
		}
	}
};

}



#endif /* SRC_CONFIG_CONFIGURATION_H_ */
