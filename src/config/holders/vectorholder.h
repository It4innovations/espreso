
#ifndef SRC_CONFIG_HOLDERS_VECTORHOLDER_H_
#define SRC_CONFIG_HOLDERS_VECTORHOLDER_H_

#include "config/configuration.h"
#include "basis/utilities/parser.h"

namespace espreso {

template <typename Ttype>
struct ECFVectorHolder: public ECFValue {
	std::vector<Ttype> &values;

	ECFVectorHolder(std::vector<Ttype> &values): values(values) {}

	std::string getValue() const
	{
		std::stringstream ss;
		if (values.size()) {
			ss << values.front();
		}
		for (size_t i = 1; i < values.size(); i++) {
			ss << ", " << values[i];
		}

		return ss.str();
	}

	bool _setValue(const std::string &values)
	{
		std::vector<std::string> parsed = Parser::split(values, ", ", true);
		for (size_t i = 0; i < parsed.size(); i++) {
			std::stringstream ssv(parsed[i]);
			Ttype tvalue;
			ssv >> tvalue;
			if (!ssv.eof() || ssv.fail()) {
				return false;
			}
			this->values.push_back(tvalue);
		}
		return true; // TODO: better format checking
	}

	virtual const void* data() const { return &values; }
};

}




#endif /* SRC_CONFIG_HOLDERS_VECTORHOLDER_H_ */
