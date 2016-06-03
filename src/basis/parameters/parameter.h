
#ifndef SRC_BASIS_PARAMETERS_PARAMETER_H_
#define SRC_BASIS_PARAMETERS_PARAMETER_H_

#include <cstdlib>

#include "../logging/logging.h"
#include "parser.h"

namespace espreso {

struct Configuration {
	std::string path;
	std::vector<std::string> nameless;
};

template <typename TValue>
struct Option {
	std::string name;
	TValue value;
	std::string description;
};

struct Envelope {
	virtual bool set(const std::string &value) =0;
	virtual std::string get() const =0;

	virtual Envelope* copy() =0;
	virtual size_t options() const
	{
		return 0;
	}

	virtual std::string optionName(size_t index) const
	{
		return "";
	}

	virtual std::string optionDesc(size_t index) const
	{
		return "";
	}

	virtual ~Envelope() {};
};

template <typename TPointer>
struct DataEnvelope: public Envelope {
	DataEnvelope(TPointer &data): _data(data) {};

	bool set(const std::string &value)
	{
		std::stringstream ss(value);
		ss >> _data;
		return ss.eof();
	}

	std::string get() const
	{
		std::stringstream ss;
		ss << _data;
		return ss.str();
	}

	Envelope* copy()
	{
		return new DataEnvelope<TPointer>(*this);
	}

private:
	TPointer &_data;
};

template <typename TPointer>
struct OptionEnvelope: public Envelope {
	OptionEnvelope(TPointer &data, std::vector<Option<TPointer> > &options)
	: _data(data), _options(options) { };

	bool set(const std::string &value)
	{
		for (size_t i = 0; i < _options.size(); i++) {
			if (StringCompare::caseInsensitiveEq(value, _options[i].name)) {
				_data = _options[i].value;
				return true;
			}
		}
		std::stringstream ss(value);
		size_t number;
		ss >> number;
		if (ss.eof() && number < _options.size()) {
			_data = _options[number].value;
			return true;
		}
		return false;
	}

	std::string get() const
	{
		for (size_t i = 0; i < _options.size(); i++) {
			if (_options[i].value == _data) {
				return _options[i].name;
			}
		}
		// this code is not reachable
		return "Unrecognized value";
	}

	Envelope* copy()
	{
		return new OptionEnvelope<TPointer>(*this);
	}

	size_t options() const
	{
		return _options.size();
	}

	std::string optionName(size_t index) const
	{
		return _options[index].name;
	}

	std::string optionDesc(size_t index) const
	{
		return _options[index].description;
	}

private:
	TPointer &_data;
	std::vector<Option<TPointer> > _options;
};

template <typename TPointer>
struct StringOptionEnvelope: public Envelope {
	StringOptionEnvelope(TPointer &data, std::vector<std::pair<std::string, std::string> > &options)
	: _data(data), _options(options) { };

	bool set(const std::string &value)
	{
		for (size_t i = 0; i < _options.size(); i++) {
			if (StringCompare::caseInsensitiveEq(value, _options[i].first)) {
				_data = i;
				return true;
			}
		}
		std::stringstream ss(value);
		size_t number;
		ss >> number;
		if (ss.eof() && number < _options.size()) {
			_data = number;
			return true;
		}
		return false;
	}

	std::string get() const
	{
		return _options[_data].first;
	}

	Envelope* copy()
	{
		return new StringOptionEnvelope<TPointer>(*this);
	}

	size_t options() const
	{
		return _options.size();
	}

	std::string optionName(size_t index) const
	{
		return _options[index].first;
	}

	std::string optionDesc(size_t index) const
	{
		return _options[index].second;
	}

private:
	TPointer &_data;
	std::vector<std::pair<std::string, std::string> > _options;
};

struct Parameter {
	enum class Help {
		INGNORE,
		WRITE
	};

	std::string name;
	Envelope *data;
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

	template <typename TParameter>
	Parameter(std::string name, TParameter &defaultValue, std::string description, Help writeToHelp = Help::INGNORE)
	: name(name), description(description), help(writeToHelp)
	{
		data = new DataEnvelope<TParameter>(defaultValue);
	}

	template <typename TParameter>
	Parameter(std::string name, TParameter &defaultValue, std::string description, std::vector<std::pair<std::string, std::string> > options, Help writeToHelp = Help::INGNORE)
	: name(name), description(description), help(writeToHelp)
	{
		data = new StringOptionEnvelope<TParameter>(defaultValue, options);
	}

	template <typename TParameter>
	Parameter(std::string name, TParameter &defaultValue, std::string description, std::vector<Option<TParameter> > options, Help writeToHelp = Help::INGNORE)
	: name(name), description(description), help(writeToHelp)
	{
		data = new OptionEnvelope<TParameter>(defaultValue, options);
	}

	Parameter(const Parameter &other)
	{
		name = other.name;
		data = other.data->copy();
		description = other.description;
		help = other.help;
	}

	Parameter& operator=(const Parameter &other)
	{
		if (this != &other) {
			name = other.name;
			data = other.data->copy();
			description = other.description;
			help = other.help;
		}
		return *this;
	}

	~Parameter()
	{
		delete data;
	}

	void set(const std::string &value)
	{
		if (!data->set(value)) {
			ESINFO(GLOBAL_ERROR) << "Parameter '" << name << "' has a wrong value '" << value << "'.";
		}
	}

	std::string get() const
	{
		return data->get();
	}
};

}


#endif /* SRC_BASIS_PARAMETERS_PARAMETER_H_ */
