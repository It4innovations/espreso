
#ifndef SRC_BASIS_PARAMETERS_ENVELOPES_STRINGOPTIONS_H_
#define SRC_BASIS_PARAMETERS_ENVELOPES_STRINGOPTIONS_H_

#include "envelope.h"

namespace espreso {

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

	void* value() const
	{
		return &_data;
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

}

#endif /* SRC_BASIS_PARAMETERS_ENVELOPES_STRINGOPTIONS_H_ */
