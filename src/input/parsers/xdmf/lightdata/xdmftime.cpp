
#include "xdmftime.h"
#include "esinfo/eslog.hpp"

#include <sstream>

using namespace espreso;

XDMFTime::XDMFTime()
: type(Type::Single), value(0)
{

}

void XDMFTime::parse(XML::Element *e)
{
	for (auto attr = e->attributes.begin(); attr != e->attributes.end(); ++attr) {
		bool _name = false, _value = false;
		if (StringCompare::caseInsensitiveEq(attr->first, "Type") || StringCompare::caseInsensitiveEq(attr->first, "TimeType")) {
			_name = true;
			_value |= set(attr->second, "Single", type, Type::Single);
			_value |= set(attr->second, "HyperSlab", type, Type::HyperSlab);
			_value |= set(attr->second, "List", type, Type::List);
			_value |= set(attr->second, "Range", type, Type::Range);
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "Value")) {
			_name = _value = true;
			std::stringstream ss(attr->second);
			ss >> value;
		}
		if (_name == false) {
			eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFTime.\n", attr->first.c_str(), attr->second.c_str());
		}
		if (_name == true && _value == false) {
			eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFTime.\n", attr->first.c_str());
		}
	}
	XDMFElement::parse(e);
}
