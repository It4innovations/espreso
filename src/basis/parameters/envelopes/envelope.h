
#ifndef SRC_BASIS_PARAMETERS_ENVELOPES_ENVELOPE_H_
#define SRC_BASIS_PARAMETERS_ENVELOPES_ENVELOPE_H_

#include "../parser.h"

namespace espreso {

template <typename TValue>
struct Option2 {
	std::string name;
	TValue value;
	std::string description;
};

struct Envelope {
	virtual bool set(const std::string &value) =0;
	virtual bool set(const std::string &attribute, const std::string &value) { return false; };
	virtual std::string get() const =0;
	virtual std::string get(const std::string &attribute) const { return ""; };

	virtual std::vector<std::string> attributes() const { return std::vector<std::string>(); };
	virtual void* value() const =0;

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

}

#endif /* SRC_BASIS_PARAMETERS_ENVELOPES_ENVELOPE_H_ */
