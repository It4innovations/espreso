
#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>

#include "configuration.h"

enum DataType {
	STRING_PARAMETER,
	INTEGER_PARAMETER,
	BOOLEAN_PARAMETER
};

class Parameter {

public:
	bool match(std::string &line) const
	{
		return line.compare(0, _name.size(), _name) == 0;
	}

	std::string name()
	{
		return _name;
	}

	std::string description()
	{
		return _description;
	}

	DataType type()
	{
		return _type;
	}

	virtual void set(std::string &line) =0;

	void error(std::string message) const
	{
		std::cerr << "The configuration file is not valid: " << message << "\n";
		exit(EXIT_FAILURE);
	}

protected:
	Parameter(DataType type, std::string name, std::string description)
		: _type(type), _delimiter("="), _name(name), _description(description) {};
	virtual ~Parameter() {};

	std::string value(std::string &line)
	{
		size_t pos = line.find(_delimiter);
		if (pos == std::string::npos) {
			error("Incorrect format of " + _name + ". Use " + _name + _delimiter + "value.");
		}
		std::string val = line.substr(pos + 1);
		return val.erase(0, val.find_first_not_of(" "));
	}

	std::string _delimiter;
	std::string _name;
	std::string _description;
	DataType _type;

};

class StringParameter : public Parameter {

public:
	StringParameter(std::string name, std::string description, std::string value)
		:Parameter(STRING_PARAMETER, name, description), _value(value) {};

	void set(std::string &line)
	{
		_value = value(line);
	}

	std::string get() { return _value; };

private:
	std::string _value;
};

class IntegerParameter : public Parameter {

public:
	IntegerParameter(std::string name, std::string description, eslocal value)
		:Parameter(INTEGER_PARAMETER, name, description), _value(0) {};

	void set(std::string &line)
	{
		std::stringstream ss(value(line));
		ss >> _value;
	}

	eslocal get() { return _value; };

private:
	eslocal _value;
};

class BooleanParameter : public Parameter {

public:
	BooleanParameter(std::string name, std::string description, eslocal value)
		:Parameter(BOOLEAN_PARAMETER, name, description), _value(0) {};

	void set(std::string &line)
	{
		size_t pos = line.find(_delimiter);
		if (pos == std::string::npos) {
			_value = true;
		} else {
			error("Boolean parameter " + _name + " should be without assignment.");
		}
	}

	bool get() { return _value; };

private:
	bool _value;
};





#endif /* PARAMETER_H_ */
