
#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "parameter.h"

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

namespace esinput {

template<class TSettings>
class Configuration {

public:
	Configuration(int argc, char** argv);
	~Configuration();

	void print() const;

	template<class ParameterType>
	const ParameterType& value(const std::string &parameter, ParameterType defaultValue) const
	{
		return _getValue(parameter, defaultValue);
	}

private:

	const eslocal& _getValue(const std::string &parameter, eslocal defaultValue) const
	{
		if (_parameters.find(parameter) != _parameters.end()) {
			return static_cast<IntegerParameter*>(_parameters.find(parameter)->second)->get();
		} else {
			std::cerr << "Incorrect parameter '" << parameter << "'.\n";
			exit(EXIT_FAILURE);
		}
	}

	const double& _getValue(const std::string &parameter, double defaultValue) const
	{
		if (_parameters.find(parameter) != _parameters.end()) {
			return static_cast<DoubleParameter*>(_parameters.find(parameter)->second)->get();
		} else {
			std::cerr << "Incorrect parameter '" << parameter << "'.\n";
			exit(EXIT_FAILURE);
		}
	}

	const std::string& _getValue(const std::string &parameter, std::string &defaultValue) const
	{
		if (_parameters.find(parameter) != _parameters.end()) {
			return static_cast<StringParameter*>(_parameters.find(parameter)->second)->get();
		} else {
			std::cerr << "Incorrect parameter '" << parameter << "'.\n";
			exit(EXIT_FAILURE);
		}
	}

	const bool& _getValue(const std::string &parameter, bool defaultValue) const
	{
		if (_parameters.find(parameter) != _parameters.end()) {
			return static_cast<BooleanParameter*>(_parameters.find(parameter)->second)->get();
		} else {
			std::cerr << "Incorrect parameter '" << parameter << "'.\n";
			exit(EXIT_FAILURE);
		}
	}

	void load(int argc, char** argv);

	std::map<std::string, Parameter*> _parameters;
};

}

#include "configuration.hpp"

#endif /* CONFIGURATION_H_ */
