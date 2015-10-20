
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

class Configuration {

public:
	Configuration(std::vector<Description> &description, int argc, char** argv);
	~Configuration();

	void print() const;

	bool isSet(const std::string &parameter) const
	{
		if (_parameters.find(parameter) != _parameters.end()) {
			return _parameters.find(parameter)->second->isSet();
		} else {
			return false;
		}
	}

	template<class ParameterType>
	const ParameterType value(const std::string &parameter, ParameterType defaultValue) const
	{
		if (_parameters.find(parameter) != _parameters.end()) {
			return _getValue(parameter, defaultValue);
		} else {
			return defaultValue;
		}
	}

private:

	eslocal _getValue(const std::string &parameter, eslocal defaultValue) const;
	double _getValue(const std::string &parameter, double defaultValue) const;
	std::string _getValue(const std::string &parameter, std::string &defaultValue) const;
	bool _getValue(const std::string &parameter, bool defaultValue) const;

	void load(int argc, char** argv);

	std::map<std::string, Parameter*> _parameters;
};

}

#endif /* CONFIGURATION_H_ */
