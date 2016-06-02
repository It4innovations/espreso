
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
	std::vector<std::pair<std::string, std::string> > options;
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

	Parameter(std::string name, int &defaultValue,         std::string description, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, long &defaultValue,        std::string description, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, size_t &defaultValue,      std::string description, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, double &defaultValue,      std::string description, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, std::string &defaultValue, std::string description, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, bool &defaultValue,        std::string description, Help writeToHelp = Help::INGNORE);


	Parameter(std::string name, int &defaultValue,         std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, long &defaultValue,        std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, size_t &defaultValue,      std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, double &defaultValue,      std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, std::string &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp = Help::INGNORE);
	Parameter(std::string name, bool &defaultValue,        std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp = Help::INGNORE);

	void set(const char* value);
	void set(const std::string &value)
	{
		set(value.c_str());
	}

	std::string get() const;
};

}


#endif /* SRC_BASIS_PARAMETERS_PARAMETER_H_ */
