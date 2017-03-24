
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

#define SUBMAPTOCONFIG(type1, type2, name, description, dName, dDescription) \
	std::map<type1, type2*> name = mapToConfiguration<type1, type2>::create(#name, description, this, dName, dDescription)

#define SUBMAPTOMAPTOCONFIG(type1, type2, type3, name, description, d1Name, d1Description, d2Name, d2Description) \
	std::map<type1, std::map<type2, type3*> > name = mapToMapToConfiguration<type1, type2, type3>::create(#name, description, this, d1Name, d1Description, d2Name, d2Description)


#define SUBMAP(type1, type2, name, description, dParameter, dValue) \
	std::map<type1, type2> name = mapToBaseType<type1, type2>::create(#name, description, this, dParameter, dValue, #type2)

#define SUBMAPTOMAP(type1, type2, type3, name, description, dParameter1, dDescription, dParameter2, dValue) \
	std::map<type1, std::map<type2, type3> > name = mapToMapToBaseType<type1, type2, type3>::create(#name, description, this, dParameter1, dDescription, dParameter2, dValue, #type3)


#define SUBMULTIMAP(type1, type2, name, description, dParameter, dValue) \
	std::multimap<type1, type2> name = multimapToBaseType<type1, type2>::create(#name, description, this, dParameter, dValue, #type2)


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
	Tvalue dummyHolder;
	std::vector<Configuration*> dummy = { &dummyHolder };

	static std::map<Tparameter, Tvalue*> create(
			const std::string &name, const std::string &description, Configuration* conf,
			const std::string &dummyName, const std::string &dummyDesc)
	{
		std::map<Tparameter, Tvalue*> configuration;
		mapToConfiguration<Tparameter, Tvalue> *subconf = new mapToConfiguration<Tparameter, Tvalue>();
		if (conf != NULL) {
			conf->subconfigurations[name] = subconf;
			conf->orderedSubconfiguration.push_back(subconf);
			conf->toDelete.push_back(subconf);
		}
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
		subconf->dummyHolder.name = dummyName;
		subconf->dummyHolder.description = dummyDesc;
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
};

template <typename Tparameter1, typename Tparameter2, typename Tvalue>
struct mapToMapToConfiguration: public Configuration {
	std::map<Tparameter1, std::map<Tparameter2, Tvalue*> > *map;
	std::map<Tparameter1, mapToConfiguration<Tparameter2, Tvalue> > submap;
	mapToConfiguration<Tparameter2, Tvalue> dummyHolder;
	std::vector<Configuration*> dummy = { &dummyHolder };

	static std::map<Tparameter1, std::map<Tparameter2, Tvalue*> > create(
			const std::string &name, const std::string &description, Configuration* conf,
			const std::string &d1Name, const std::string &d1Desc,
			const std::string &d2Name, const std::string &d2Desc)
	{
		std::map<Tparameter1, std::map<Tparameter2, Tvalue*> > configuration;
		mapToMapToConfiguration<Tparameter1, Tparameter2, Tvalue> *subconf = new mapToMapToConfiguration<Tparameter1, Tparameter2, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
		subconf->dummyHolder.name = d1Name;
		subconf->dummyHolder.description = d1Desc;
		subconf->dummyHolder.dummyHolder.name = d2Name;
		subconf->dummyHolder.dummyHolder.description = d2Desc;
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
};

template <typename Tparameter, typename Tvalue>
struct mapToBaseType: public Configuration {
	std::map<Tparameter, Tvalue> *map;
	std::string dParameter;
	Tvalue dummyHolder;
	std::vector<ParameterBase*> dummy;

	static std::map<Tparameter, Tvalue> create(
			const std::string &name, const std::string &description, Configuration* conf,
			const std::string &dParameter, const Tvalue &dValue, const std::string &type)
	{
		std::map<Tparameter, Tvalue> configuration;
		mapToBaseType<Tparameter, Tvalue> *subconf = new mapToBaseType<Tparameter, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
		subconf->dummyHolder = dValue;
		subconf->dummy.push_back(new ValueHolder<Tvalue>(dParameter, "Accepts list of parameters of the following type: " + dParameter, subconf->dummyHolder, subconf->dummyHolder, type));
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

template <typename Tparameter, typename Tvalue>
struct multimapToBaseType: public Configuration {
	std::multimap<Tparameter, Tvalue> *map;
	std::string dParameter;
	Tvalue dummyHolder;
	std::vector<ParameterBase*> dummy;

	static std::multimap<Tparameter, Tvalue> create(
			const std::string &name, const std::string &description, Configuration* conf,
			const std::string &dParameter, const Tvalue &dValue, const std::string &type)
	{
		std::multimap<Tparameter, Tvalue> configuration;
		multimapToBaseType<Tparameter, Tvalue> *subconf = new multimapToBaseType<Tparameter, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
		subconf->dummyHolder = dValue;
		subconf->dummy.push_back(new ValueHolder<Tvalue>(dParameter, "Accepts list of (possible same) parameters of the following type: " + dParameter, subconf->dummyHolder, subconf->dummyHolder, type));
		return configuration;
	}

	bool set(const std::string &parameter, const std::string &value)
	{
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
		auto it = map->insert(std::make_pair(pholder.value, vholder.value));
		orderedParameters.push_back(new ValueHolder<Tvalue>(parameter, "Parameter value", it->second, it->second, ""));
		return true;
	}

	virtual const std::vector<ParameterBase*>& storeParameters() const
	{
		if (orderedParameters.size()) {
			return orderedParameters;
		} else {
			return dummy;
		}
	}

	~multimapToBaseType()
	{
		if (!copy) {
			std::for_each(dummy.begin(), dummy.end(), [] (ParameterBase * p) { delete p; });
			std::for_each(orderedParameters.begin(), orderedParameters.end(), [] (ParameterBase * p) { delete p; });
		}
	}
};

template <typename Tparameter1, typename Tparameter2, typename Tvalue>
struct mapToMapToBaseType: public Configuration {
	std::map<Tparameter1, std::map<Tparameter2, Tvalue> > *map;
	std::map<Tparameter1, mapToBaseType<Tparameter2, Tvalue> > submap;
	mapToBaseType<Tparameter2, Tvalue> dummyHolder;
	std::vector<Configuration*> dummy = { &dummyHolder };

	static std::map<Tparameter1, std::map<Tparameter2, Tvalue> > create(
			const std::string &name, const std::string &description, Configuration* conf,
			const std::string &dParameter1, const std::string &dDescription,
			const std::string &dParameter2, const std::string &dValue, const std::string &type)
	{
		std::map<Tparameter1, std::map<Tparameter2, Tvalue> > configuration;
		mapToMapToBaseType<Tparameter1, Tparameter2, Tvalue> *subconf = new mapToMapToBaseType<Tparameter1, Tparameter2, Tvalue>();
		conf->subconfigurations[name] = subconf;
		conf->orderedSubconfiguration.push_back(subconf);
		conf->toDelete.push_back(subconf);
		subconf->map = &configuration;
		subconf->name = name;
		subconf->description = description;
		subconf->dummyHolder.name = dParameter1;
		subconf->dummyHolder.description = dDescription;
		subconf->dummyHolder.dummyHolder = dValue;
		subconf->dummyHolder.dummy.push_back(new ValueHolder<Tvalue>(dParameter2, "Accepts list of parameters of the following type: " + dParameter2, subconf->dummyHolder.dummyHolder, subconf->dummyHolder.dummyHolder, type));
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
