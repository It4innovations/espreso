
#ifndef SRC_CONFIG_CONFIGURATION_H_
#define SRC_CONFIG_CONFIGURATION_H_

#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>

#include "esbasis.h"
#include "reader.h"

#define MAX_VECTOR_CONFIGURATION_SIZE 10

#define OPTIONS(...) __VA_ARGS__
#define OPTION(type, name, description, value, options)         type name = DataHolder::create<type>(#name, description, name, value, options, this)
#define PARAMETER(type, name, description, value)               type name = DataHolder::create<type>(#name, description, name, value, #type, this)

#define SUBVECTOR(type, name, description, dparameter, dvalue)           ConfigurationVector<type> name = ConfigurationVector<type>::create(#name, description, dparameter, dvalue, this)
#define SUBMAP(ptype, vtype, name, description, dparameter, dvalue) ConfigurationMap<ptype, vtype> name = ConfigurationMap<ptype, vtype>::create(#name, description, #vtype, dparameter, dvalue, this)
#define SUBCONFIG(type, name, description)                                                    type name = Configuration::create<type>(#name, description, this)

namespace espreso {

struct ParameterBase {
	std::string name;
	std::string description;
	std::string allowedValue;

	ParameterBase(const std::string &name, const std::string &description, const std::string &allowedValue)
	: name(name), description(description), allowedValue(allowedValue) {}

	virtual bool set(const std::string &value) =0;
	virtual std::string get() const =0;

	virtual ~ParameterBase() {};
};

template <typename Tvalue>
struct Option {
	std::string name;
	std::string description;
	Tvalue value;

	Option(const std::string &name, Tvalue value, const std::string &description): name(name), description(description), value(value) {}
};



template <typename Ttype>
struct ValueHolder: public ParameterBase {
	Ttype &value;

	ValueHolder(const std::string &name, const std::string &description, Ttype &value, Ttype defaultValue, const std::string &type)
	: ParameterBase(name, description, type), value(value) { value = defaultValue; }

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
struct ValueHolder<std::string>: public ParameterBase {
	std::string &value;

	ValueHolder(const std::string &name, const std::string &description, std::string &value, const std::string &defaultValue, const std::string &type)
	: ParameterBase(name, description, "*"), value(value) { };

	bool set(const std::string &value)
	{
		this->value = value;
		return true;
	}

	std::string get() const
	{
		return value;
	}
};

template <typename Ttype>
struct OptionsHolder: public ParameterBase {
	Ttype &value;
	std::vector<Option<Ttype> > options;

	OptionsHolder(const std::string &name, const std::string &description, Ttype &value, Ttype defaultValue, std::vector<Option<Ttype> > options)
	: ParameterBase(name, description, ""), value(value), options(options)
	{
		if (options.size()) {
			std::for_each(options.begin(), options.end() - 1, [&] (Option<Ttype> &option) { allowedValue += option.name + ", "; });
			allowedValue += options.back().name;
		}
		value = defaultValue;
	}

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

struct Configuration {
	std::map<std::string, ParameterBase*, StringCompare> parameters;
	std::map<std::string, Configuration*, StringCompare> subconfigurations;
	std::vector<ParameterBase*> orderedParameters;
	std::vector<Configuration*> orderedSubconfiguration;
	std::string name;
	std::string description;

	static void read(const std::string &file) { Reader::read(file); }
	static void read(int *argc, char ***argv) { Reader::read(argc, argv); }
	static void print() { Reader::print(); }
	static void store() { Reader::store(); }

	template <typename Ttype>
	static Ttype create(const std::string &name, const std::string &description, Configuration* conf)
	{
		Ttype configuration;
		conf->subconfigurations[name] = &configuration;
		conf->orderedSubconfiguration.push_back(&configuration);
		configuration.name = name;
		configuration.description = description;
		return configuration;
	}

	virtual Configuration& operator[](const std::string &subconfiguration)
	{
		if (subconfigurations.find(subconfiguration) != subconfigurations.end()) {
			return *subconfigurations.find(subconfiguration)->second;
		} else {
			ESINFO(GLOBAL_ERROR) << "Unknown sub-configuration '" << subconfiguration << "'";
			exit(EXIT_FAILURE);
		}
	}

	virtual bool set(const std::string &parameter, const std::string &value)
	{
		if (parameters.find(parameter) != parameters.end()) {
			return parameters.find(parameter)->second->set(value);
		} else {
			ESINFO(GLOBAL_ERROR) << "Unknown parameter '" << parameter << "'";
			return false;
		}
	}

	virtual const std::vector<Configuration*>& storeConfigurations() const { return orderedSubconfiguration; }
	virtual const std::vector<ParameterBase*>& storeParameters() const { return orderedParameters; }

	virtual ~Configuration()
	{
		for (auto it = parameters.begin(); it != parameters.end(); ++it) {
			delete it->second;
		}
	}
};

template <typename Ttype>
struct ConfigurationVector: public Configuration {
	std::vector<Configuration*> configurations;
	std::vector<Configuration*> dummy;

	static ConfigurationVector<Ttype> create(const std::string &name, const std::string &description, const std::string &dParameter, const std::string &dValue, Configuration* conf)
	{
		ConfigurationVector<Ttype> configuration;
		conf->subconfigurations[name] = &configuration;
		conf->orderedSubconfiguration.push_back(&configuration);
		configuration.name = name;
		configuration.description = description;
		configuration.dummy.push_back(new Ttype{});
		configuration.dummy.back()->name = dParameter;
		configuration.dummy.back()->description = dValue;
		return configuration;
	}

	Configuration& operator[](const std::string &subconfiguration)
	{
		std::stringstream ss(subconfiguration);
		size_t index;
		ss >> index;
		if (!ss.eof() || ss.fail()) {
			ESINFO(GLOBAL_ERROR) << "Invalid vector index '" << subconfiguration << "'";
		}
		if (index > configurations.size()) {
			configurations.resize(index, NULL);
			configurations[index - 1] = new Ttype{};
			subconfigurations[std::to_string(index - 1)] = configurations[index - 1];
			orderedSubconfiguration.push_back(configurations[index - 1]);
			configurations[index - 1]->name = std::to_string(index);
		}
		return *configurations[index - 1];
	}

	virtual const std::vector<Configuration*>& storeConfigurations() const
	{
		if (orderedSubconfiguration.size()) {
			return orderedSubconfiguration;
		} else {
			return dummy;
		}
	}

	~ConfigurationVector()
	{
		std::for_each(configurations.begin(), configurations.end(), [] (Configuration * c) { delete c; });
		std::for_each(dummy.begin(), dummy.end(), [] (Configuration * c) { delete c; });
	}
};

template <typename Tparameter, typename Tvalue>
struct ConfigurationMap: public Configuration {
	std::map<Tparameter, Tvalue> values;
	Tvalue dummyValue;
	std::vector<ParameterBase*> dummy;
	std::string type;

	static ConfigurationMap<Tparameter, Tvalue> create(const std::string &name, const std::string &description, const std::string &type, const std::string &dParameter, Tvalue dValue, Configuration* conf)
	{
		ConfigurationMap<Tparameter, Tvalue> configuration;
		conf->subconfigurations[name] = &configuration;
		conf->orderedSubconfiguration.push_back(&configuration);
		configuration.name = name;
		configuration.description = description;
		configuration.dummy.push_back(new ValueHolder<Tvalue>("# " + dParameter, "List of values.", configuration.dummyValue, dValue, type));
		configuration.dummyValue = dValue;
		return configuration;
	}

	bool set(const std::string &parameter, const std::string &value)
	{
		if (parameters.find(parameter) != parameters.end()) {
			return parameters.find(parameter)->second->set(value);
		} else {
			Tparameter par;
			ValueHolder<Tparameter> pholder(parameter, "", par, par, type);
			if (!pholder.set(parameter)) {
				return false;
			}
			Tparameter val;
			ValueHolder<Tparameter> vholder(value, "", val, val, type);
			if (!vholder.set(value)) {
				return false;
			}
			values[pholder.value] = vholder.value;
			parameters[parameter] = new ValueHolder<Tvalue>(parameter, "Parameter value", values[pholder.value], values[pholder.value], type);
			orderedParameters.push_back(parameters[parameter]);
			return true;
		}
	}

	virtual const std::vector<ParameterBase*>& storeParameters() const
	{
		if (orderedParameters.size()) {
			return orderedParameters;
		} else {
			return dummy;
		}
	}

	~ConfigurationMap()
	{
		std::for_each(dummy.begin(), dummy.end(), [] (ParameterBase * p) { delete p; });
	}
};

struct DataHolder {
	template <typename Ttype>
	static Ttype create(const std::string &name, const std::string &description, Ttype &value, Ttype defaultValue, const std::string &type, Configuration* configuration)
	{
		ParameterBase *parameter = new ValueHolder<Ttype>(name, description, value, defaultValue, type);
		configuration->parameters[name] = parameter;
		configuration->orderedParameters.push_back(parameter);
		return defaultValue;
	}

	template <typename Ttype>
	static Ttype create(const std::string &name, const std::string &description, Ttype &value, Ttype defaultValue, std::vector<Option<Ttype> > options, Configuration* configuration)
	{
		ParameterBase *parameter = new OptionsHolder<Ttype>(name, description, value, defaultValue, options);
		configuration->parameters[name] = parameter;
		configuration->orderedParameters.push_back(parameter);
		return defaultValue;
	}
};

}



#endif /* SRC_CONFIG_CONFIGURATION_H_ */
