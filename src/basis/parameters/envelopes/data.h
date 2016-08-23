
#ifndef SRC_BASIS_PARAMETERS_ENVELOPES_DATA_H_
#define SRC_BASIS_PARAMETERS_ENVELOPES_DATA_H_

#include "envelope.h"

namespace espreso {

template <typename TPointer>
struct DataEnvelope: public Envelope {
	DataEnvelope(TPointer &data): _data(data) {};

	bool set(const std::string &value)
	{
		std::stringstream ss(value);
		ss >> _data;
		return ss.eof() && !ss.fail();
	}

	std::string get() const
	{
		std::stringstream ss;
		ss << _data;
		return ss.str();
	}

	void* value() const
	{
		return &_data;
	}

	Envelope* copy()
	{
		return new DataEnvelope<TPointer>(*this);
	}

private:
	TPointer &_data;
};

template <>
struct DataEnvelope<bool>: public Envelope {
	DataEnvelope(bool &data): _data(data) {};

	bool set(const std::string &value)
	{
		if (value.size() == 0) {
			_data = true;
			return true;
		} else {
			std::stringstream ss(value);
			ss >> _data;
			return ss.eof();
		}
	}

	std::string get() const
	{
		std::stringstream ss;
		ss << _data;
		return ss.str();
	}

	void* value() const
	{
		return &_data;
	}

	Envelope* copy()
	{
		return new DataEnvelope<bool>(*this);
	}

private:
	bool &_data;
};

template <>
struct DataEnvelope<std::string>: public Envelope {
	DataEnvelope(std::string &data): _data(data) {};

	bool set(const std::string &value)
	{
		_data = value;
		return true;
	}

	std::string get() const
	{
		return _data;
	}

	void* value() const
	{
		return &_data;
	}

	Envelope* copy()
	{
		return new DataEnvelope<std::string>(*this);
	}

private:
	std::string &_data;
};

}

#endif /* SRC_BASIS_PARAMETERS_ENVELOPES_DATA_H_ */
