
#ifndef SRC_BASIS_PARAMETERS_PARAMETER_H_
#define SRC_BASIS_PARAMETERS_PARAMETER_H_

#include "../logging/logging.h"
#include "parser.h"

namespace espreso {

struct Configuration {
	std::string path;
	std::vector<std::string> nameless;
};

struct Parameter {
	enum class Help {
		INGNORE,
		WRITE
	};

	enum class DataType {
		STRING,
		INTEGER,
		LONG,
		SIZE,
		DOUBLE,
		BOOLEAN
	};

	DataType type;
	std::string name;
	void* value;
	std::string description;
	std::vector<std::string> options;
	Help help;

	bool operator<(const Parameter &other)
	{
		return StringCompare::caseInsensitive(name, other.name);
	}

	bool operator<(Parameter &other) const
	{
		return StringCompare::caseInsensitive(name, other.name);
	}

	bool operator<(Parameter &other)
	{
		return StringCompare::caseInsensitive(name, other.name);
	}

	bool operator<(const std::string &str) const
	{
		return StringCompare::caseInsensitive(name, str);
	}

	Parameter(std::string name, int &defaultValue, std::string description, Help writeToHelp = Help::INGNORE)
	: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

	Parameter(std::string name, long &defaultValue, std::string description, Help writeToHelp = Help::INGNORE)
	: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

	Parameter(std::string name, size_t &defaultValue, std::string description, Help writeToHelp = Help::INGNORE)
	: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

	Parameter(std::string name, double &defaultValue, std::string description, Help writeToHelp = Help::INGNORE)
	: type(DataType::DOUBLE), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

	Parameter(std::string name, std::string &defaultValue, std::string description, Help writeToHelp = Help::INGNORE)
	: type(DataType::STRING), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

	Parameter(std::string name, bool &defaultValue, std::string description, Help writeToHelp = Help::INGNORE)
	: type(DataType::BOOLEAN), name(name), value(&defaultValue), description(description), help(writeToHelp) { };


	Parameter(std::string name, int &defaultValue, std::string description, std::vector<std::string> options, Help writeToHelp = Help::INGNORE)
	: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

	Parameter(std::string name, long &defaultValue, std::string description, std::vector<std::string> options, Help writeToHelp = Help::INGNORE)
	: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

	Parameter(std::string name, size_t &defaultValue, std::string description, std::vector<std::string> options, Help writeToHelp = Help::INGNORE)
	: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

	Parameter(std::string name, double &defaultValue, std::string description, std::vector<std::string> options, Help writeToHelp = Help::INGNORE)
	: type(DataType::DOUBLE), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

	Parameter(std::string name, std::string &defaultValue, std::string description, std::vector<std::string> options, Help writeToHelp = Help::INGNORE)
	: type(DataType::STRING), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

	Parameter(std::string name, bool &defaultValue, std::string description, std::vector<std::string> options, Help writeToHelp = Help::INGNORE)
	: type(DataType::BOOLEAN), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

	void set(const std::string &value)
	{
		set(value.c_str());
	}

	void set(const char* value)
	{
		switch (type) {
		case DataType::INTEGER:
			*(static_cast<int*>(this->value)) = std::stoi(value);
			break;
		case DataType::LONG:
			*(static_cast<long*>(this->value)) = std::stol(value);
			break;
		case DataType::SIZE:
			*(static_cast<size_t*>(this->value)) = std::stoull(value);
			break;
		case DataType::DOUBLE:
			*(static_cast<double*>(this->value)) = std::stod(value);
			break;
		case DataType::STRING:
			*(static_cast<std::string*>(this->value)) = value;
			break;
		case DataType::BOOLEAN:
			*(static_cast<bool*>(this->value)) = (value[0] == '0' && value[1] == 0) ? false : true;
			break;
		}
	}

	std::string get() const
	{
		std::stringstream ss;
		switch (type) {
		case DataType::INTEGER:
			ss << *(static_cast<int*>(this->value));
			return ss.str();
		case DataType::LONG:
			ss << *(static_cast<long*>(this->value));
			return ss.str();
		case DataType::SIZE:
			ss << *(static_cast<size_t*>(this->value));
			return ss.str();
		case DataType::DOUBLE:
			ss << *(static_cast<double*>(this->value));
			return ss.str();
		case DataType::STRING:
			ss << *(static_cast<std::string*>(this->value));
			return ss.str();
		case DataType::BOOLEAN:
			if (*(static_cast<bool*>(this->value))) {
				return "TRUE";
			} else {
				return "FALSE";
			}
		}
		return "";
	}
};

}


#endif /* SRC_BASIS_PARAMETERS_PARAMETER_H_ */
