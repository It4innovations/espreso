
#ifndef SRC_CONFIGURATION_CONFIGURATION_HPP_
#define SRC_CONFIGURATION_CONFIGURATION_HPP_

#include <sstream>
#include <algorithm>

#include "configuration.h"
#include "../basis/logging/logging.h"

#define OPTIONS(...) __VA_ARGS__
#define OPTION(type, name, description, value, options)         type name = ParameterHolder::create<type>(#name, description, name, value, options, this)
#define PARAMETER(type, name, description, value)               type name = ParameterHolder::create<type>(#name, description, name, value, #type, this)

#define SUBCONFIG(type, name, description) \
	type name = Configuration::create<type>(#name, description, this)

#define SUBMAP(type1, type2, name, description, dparameter, dvalue) \
	std::map<type1, type2> name = mapToBaseType<type1, type2>::create(#name, description, this)

#define SUBMAPTOMAP(type1, type2, type3, name, description) \
	std::map<type1, std::map<type2, type3> > name = mapToMapToBaseType<type1, type2, type3>::create(#name, description, this)

#define SUBMAPTOCONFIG(type1, type2, name, description) \
	std::map<type1, type2*> name = mapToConfiguration<type1, type2>::create(#name, description, this)

#define SUBMAPTOMAPTOCONFIG(type1, type2, type3, name, description) \
	std::map<type1, std::map<type2, type3*> > name = mapToMapToConfiguration<type1, type2, type3>::create(#name, description, this)

namespace espreso {

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

	size_t index() const { return (size_t)value; }
};

template <typename Tparameter, typename Tvalue>
struct mapToConfiguration: public Configuration {
	std::map<Tparameter, Tvalue*> *map;
	std::vector<Configuration*> dummy;

	static std::map<Tparameter, Tvalue*> create(const std::string &name, const std::string &description, Configuration* conf)
	{
		std::map<Tparameter, Tvalue*> configuration;
		mapToConfiguration<Tparameter, Tvalue> *subconf = new mapToConfiguration<Tparameter, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
		subconf->dummy.push_back(new mapToConfiguration<Tparameter, Tvalue>());
		subconf->dummy.back()->name = "DEFAULT";
		subconf->dummy.back()->description = "VALUE"; // TODO: improve
		return configuration;
	}

	Configuration& operator[](const std::string &parameter)
	{
		Tparameter param;
		ValueHolder<Tparameter> pholder(parameter, "", param, param, "");
		if (!pholder.set(parameter)) {
			ESINFO(GLOBAL_ERROR) << "Invalid object type of object '" << name << "'";
		}

		if (map->find(param) == map->end()) {
			Tvalue *value = new Tvalue{};
			(*map)[param] = value;
			toDelete.push_back(value);
			subconfigurations[parameter] = value;
			orderedSubconfiguration.push_back(value);
			subconfigurations[parameter]->name = parameter;
		}
		return *map->find(param)->second;
	}

	virtual const std::vector<Configuration*>& storeConfigurations() const
	{
		if (orderedSubconfiguration.size()) {
			return orderedSubconfiguration;
		} else {
			return dummy;
		}
	}

	~mapToConfiguration()
	{
		if (!copy) {
			std::for_each(dummy.begin(), dummy.end(), [] (Configuration * c) { delete c; });
		}
	}
};

template <typename Tparameter1, typename Tparameter2, typename Tvalue>
struct mapToMapToConfiguration: public Configuration {
	std::map<Tparameter1, std::map<Tparameter2, Tvalue*> > *map;
	std::map<Tparameter1, mapToConfiguration<Tparameter2, Tvalue> > submap;
	std::vector<Configuration*> dummy;

	static std::map<Tparameter1, std::map<Tparameter2, Tvalue*> > create(const std::string &name, const std::string &description, Configuration* conf)
	{
		std::map<Tparameter1, std::map<Tparameter2, Tvalue*> > configuration;
		mapToMapToConfiguration<Tparameter1, Tparameter2, Tvalue> *subconf = new mapToMapToConfiguration<Tparameter1, Tparameter2, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
//		subconf->dummy.push_back(new mapToConfiguration<Tparameter1, Tparameter2, Tvalue>());
//		subconf->dummy.back()->name = "DEFAULT";
//		subconf->dummy.back()->description = "VALUE"; // TODO: improve
		return configuration;
	}

	Configuration& operator[](const std::string &parameter)
	{
		Tparameter1 param;
		ValueHolder<Tparameter1> pholder(parameter, "", param, param, "");
		if (!pholder.set(parameter)) {
			ESINFO(GLOBAL_ERROR) << "Invalid object type of object '" << name << "'";
		}

		if (submap.find(param) == submap.end()) {
			(*map)[param] = std::map<Tparameter2, Tvalue*>();
			auto &value = (*map)[param];
			submap[param].map = &value;
		}
		return submap.find(param)->second;
	}

	virtual const std::vector<Configuration*>& storeConfigurations() const
	{
		if (orderedSubconfiguration.size()) {
			return orderedSubconfiguration;
		} else {
			return dummy;
		}
	}

	~mapToMapToConfiguration()
	{
		if (!copy) {
			std::for_each(dummy.begin(), dummy.end(), [] (Configuration * c) { delete c; });
		}
	}
};

template <typename Tparameter, typename Tvalue>
struct mapToBaseType: public Configuration {
	std::map<Tparameter, Tvalue> *map;
	std::vector<ParameterBase*> dummy;

	static std::map<Tparameter, Tvalue> create(const std::string &name, const std::string &description, Configuration* conf)
	{
		std::map<Tparameter, Tvalue> configuration;
		mapToBaseType<Tparameter, Tvalue> *subconf = new mapToBaseType<Tparameter, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
//		subconf->dummy.push_back(new mapToConfiguration<Tparameter, Tvalue>());
//		subconf->dummy.back()->name = "DEFAULT";
//		subconf->dummy.back()->description = "VALUE"; // TODO: improve
		return configuration;
	}

	bool set(const std::string &parameter, const std::string &value)
	{
		if (parameters.find(parameter) != parameters.end()) {
			return parameters.find(parameter)->second->set(value);
		} else {
			Tparameter par;
			ValueHolder<Tparameter> pholder(parameter, "", par, par, "");
			if (!pholder.set(parameter)) {
				return false;
			}
			Tvalue val;
			ValueHolder<Tvalue> vholder(value, "", val, val, "");
			if (!vholder.set(value)) {
				return false;
			}
			(*map)[pholder.value] = vholder.value;
			parameters[parameter] = new ValueHolder<Tvalue>(parameter, "Parameter value", (*map)[pholder.value], (*map)[pholder.value], "");
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

	~mapToBaseType()
	{
		if (!copy) {
			std::for_each(dummy.begin(), dummy.end(), [] (ParameterBase * p) { delete p; });
		}
	}
};

template <typename Tparameter1, typename Tparameter2, typename Tvalue>
struct mapToMapToBaseType: public Configuration {
	std::map<Tparameter1, std::map<Tparameter2, Tvalue> > *map;
	std::map<Tparameter1, mapToBaseType<Tparameter2, Tvalue> > submap;
	std::vector<Configuration*> dummy;

	static std::map<Tparameter1, std::map<Tparameter2, Tvalue> > create(const std::string &name, const std::string &description, Configuration* conf)
	{
		std::map<Tparameter1, std::map<Tparameter2, Tvalue> > configuration;
		mapToMapToBaseType<Tparameter1, Tparameter2, Tvalue> *subconf = new mapToMapToBaseType<Tparameter1, Tparameter2, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
//		subconf->dummy.push_back(new mapToConfiguration<Tparameter1, Tparameter2, Tvalue>());
//		subconf->dummy.back()->name = "DEFAULT";
//		subconf->dummy.back()->description = "VALUE"; // TODO: improve
		return configuration;
	}

	Configuration& operator[](const std::string &parameter)
	{
		Tparameter1 param;
		ValueHolder<Tparameter1> pholder(parameter, "", param, param, "");
		if (!pholder.set(parameter)) {
			ESINFO(GLOBAL_ERROR) << "Invalid object type of object '" << name << "'";
		}

		if (submap.find(param) == submap.end()) {
			auto &value = (*map)[param];
			submap[param].map = &value;
		}
		return submap.find(param)->second;
	}

	virtual const std::vector<Configuration*>& storeConfigurations() const
	{
		if (orderedSubconfiguration.size()) {
			return orderedSubconfiguration;
		} else {
			return dummy;
		}
	}

	~mapToMapToBaseType()
	{
		if (!copy) {
			std::for_each(dummy.begin(), dummy.end(), [] (Configuration * c) { delete c; });
		}
	}
};


struct ParameterHolder {
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


#endif /* SRC_CONFIGURATION_CONFIGURATION_HPP_ */
