
#include "parameter.h"

using namespace espreso;

Parameter::Parameter(std::string name, int &defaultValue, std::string description, Help writeToHelp)
: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

Parameter::Parameter(std::string name, long &defaultValue, std::string description, Help writeToHelp)
: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

Parameter::Parameter(std::string name, size_t &defaultValue, std::string description, Help writeToHelp)
: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

Parameter::Parameter(std::string name, double &defaultValue, std::string description, Help writeToHelp)
: type(DataType::DOUBLE), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

Parameter::Parameter(std::string name, std::string &defaultValue, std::string description, Help writeToHelp)
: type(DataType::STRING), name(name), value(&defaultValue), description(description), help(writeToHelp) { };

Parameter::Parameter(std::string name, bool &defaultValue, std::string description, Help writeToHelp)
: type(DataType::BOOLEAN), name(name), value(&defaultValue), description(description), help(writeToHelp) { };


Parameter::Parameter(std::string name, int &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp)
: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

Parameter::Parameter(std::string name, long &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp)
: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

Parameter::Parameter(std::string name, size_t &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp)
: type(DataType::INTEGER), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

Parameter::Parameter(std::string name, double &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp)
: type(DataType::DOUBLE), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

Parameter::Parameter(std::string name, std::string &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp)
: type(DataType::STRING), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };

Parameter::Parameter(std::string name, bool &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp)
: type(DataType::BOOLEAN), name(name), value(&defaultValue), description(description), options(options), help(writeToHelp) { };


void Parameter::set(const char* value)
{
	size_t option = 0;
	for (size_t i = 0; i < options.size(); i++) {
		if (StringCompare::caseInsensitiveEq(options[i].first, value)) {
			option = i;
			break;
		}
	}

	switch (type) {
	case DataType::INTEGER:
		*(static_cast<int*>(this->value)) = option < options.size() ? option : std::stoi(value);
		break;
	case DataType::LONG:
		*(static_cast<long*>(this->value)) = option < options.size() ? option : std::stol(value);;
		break;
	case DataType::SIZE:
		*(static_cast<size_t*>(this->value)) = option < options.size() ? option : std::stoull(value);
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

std::string Parameter::get() const
{
	std::stringstream ss;
	switch (type) {
	case DataType::INTEGER:
		if (*(static_cast<int*>(this->value)) < options.size()) {
			ss << options[*(static_cast<int*>(this->value))].first;
		} else {
			ss << *(static_cast<int*>(this->value));
		}
		return ss.str();
	case DataType::LONG:
		if (*(static_cast<long*>(this->value)) < options.size()) {
			ss << options[*(static_cast<long*>(this->value))].first;
		} else {
			ss << *(static_cast<long*>(this->value));
		}
		return ss.str();
	case DataType::SIZE:
		if (*(static_cast<size_t*>(this->value)) < options.size()) {
			ss << options[*(static_cast<size_t*>(this->value))].first;
		} else {
			ss << *(static_cast<size_t*>(this->value));
		}
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



