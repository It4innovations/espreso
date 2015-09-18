
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
	bool match(const std::string &line) const
	{
		std::string param = line.substr(0, line.find(" "));
		return param.size() == _name.size() && param.compare(0, _name.size(), _name) == 0;
	}

	const std::string& name() const
	{
		return _name;
	}

	const std::string& description() const
	{
		return _description;
	}

	DataType type() const
	{
		return _type;
	}

	bool isSet() const
	{
		return _set;
	}

	void reset(bool value)
	{
		_set = value;
	}

	virtual void set(const std::string &line) =0;
	virtual Parameter* copy() =0;

	void error(std::string message) const
	{
		std::cerr << "The configuration file is not valid: " << message << "\n";
		exit(EXIT_FAILURE);
	}

	virtual ~Parameter() {};

protected:
	Parameter(DataType type, std::string name, std::string description)
		: _type(type), _delimiter("="), _name(name), _description(description), _set(false) {};

	std::string value(const std::string &line)
	{
		size_t pos = line.find(_delimiter);
		if (pos == std::string::npos) {
			error("Incorrect format of " + _name + ". Use " + _name + _delimiter + "value.");
		}
		std::string val = line.substr(pos + 1);
		val.erase(0, val.find_first_not_of(" "));
		if (val.find_last_of(" ") != std::string::npos) {
			return val.erase(val.find_last_of(" "));
		} else {
			return val;
		}
	}

	std::string _delimiter;
	std::string _name;
	std::string _description;
	DataType _type;
	bool _set;

};

class StringParameter : public Parameter {

public:
	StringParameter(std::string name, std::string description, std::string value)
		:Parameter(STRING_PARAMETER, name, description), _value(value) {};

	void set(const std::string &line)
	{
		_value = value(line);
		if (!_value.size()) {
			error("Empty parameter " + _name + ".\n");
		}
		_set = true;
	}

	const std::string& get() const { return _value; };

	Parameter* copy()
	{
		return new StringParameter(*this);
	}

private:
	std::string _value;
};

class IntegerParameter : public Parameter {

public:
	IntegerParameter(std::string name, std::string description, eslocal value)
		:Parameter(INTEGER_PARAMETER, name, description), _value(value) {};

	void set(const std::string &line)
	{
		std::stringstream ss(value(line));
		ss >> _value;
		_set = true;
	}

	const eslocal& get() const { return _value; };

	Parameter* copy()
	{
		return new IntegerParameter(*this);
	}

private:
	eslocal _value;
};

class BooleanParameter : public Parameter {

public:
	BooleanParameter(std::string name, std::string description, eslocal value)
		:Parameter(BOOLEAN_PARAMETER, name, description), _value(value) {};

	void set(const std::string &line)
	{
		size_t pos = line.find(_delimiter);
		if (pos == std::string::npos) {
			_value = true;
			_set = true;
		} else {
			error("Boolean parameter " + _name + " should be without assignment.");
		}
	}

	const bool& get() const { return _value; };

	Parameter* copy()
	{
		return new BooleanParameter(*this);
	}

private:
	bool _value;
};





#endif /* PARAMETER_H_ */
