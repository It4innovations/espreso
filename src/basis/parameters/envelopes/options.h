
#ifndef SRC_BASIS_PARAMETERS_ENVELOPES_OPTIONS_H_
#define SRC_BASIS_PARAMETERS_ENVELOPES_OPTIONS_H_

#include "envelope.h"

namespace espreso {

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

}

#endif /* SRC_BASIS_PARAMETERS_ENVELOPES_OPTIONS_H_ */
