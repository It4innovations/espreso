
#ifndef SRC_BASIS_PARAMETERS_ENVELOPES_MAP_H_
#define SRC_BASIS_PARAMETERS_ENVELOPES_MAP_H_

#include "envelope.h"

namespace espreso {

template <typename TPointer>
struct MapEnvelope: public Envelope {
	MapEnvelope(std::map<std::string, TPointer> &data): _data(data) {};

	bool set(const std::string &value)
	{
		return false;
	}

	std::string get() const
	{
		return "{}";
	}

	bool set(const std::string &attribute, const std::string &value)
	{
		std::stringstream ss(value);
		ss >> _data[attribute];
		return ss.eof() && !ss.fail();
	}

	std::string get(const std::string &attribute) const
	{
		std::stringstream ss;
		ss << _data[attribute];
		return ss.str();
	}

	std::vector<std::string> attributes() const
	{
		std::vector<std::string> attributes;
		for (auto it = _data.begin(); it != _data.end(); ++it) {
			attributes.push_back(it->first);
		}
		return attributes;
	}

	void* value() const
	{
		return &_data;
	}

	Envelope* copy()
	{
		return new MapEnvelope<TPointer>(*this);
	}

private:
	std::map<std::string, TPointer> &_data;
};

}




#endif /* SRC_BASIS_PARAMETERS_ENVELOPES_MAP_H_ */
