
#include "xdmfattribute.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

XDMFAttribute::XDMFAttribute()
: type(Type::Scalar), center(Center::Node), itemtype(), efaminly(), edegree(), ecell()
{

}

void XDMFAttribute::parse(XML::Element *e)
{
	for (auto attr = e->attributes.begin(); attr != e->attributes.end(); ++attr) {
		bool _name = false, _value = false;
		if (StringCompare::caseInsensitiveEq(attr->first, "Name")) {
			_name = _value = true;
			name = attr->second;
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "Reference")) {
			_name = _value = true;
			reference = attr->second;
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "Type") || StringCompare::caseInsensitiveEq(attr->first, "AttributeType")) {
			_name = true;
			_value |= set(attr->second, "Scalar", type, Type::Scalar);
			_value |= set(attr->second, "Vector", type, Type::Vector);
			_value |= set(attr->second, "Tensor", type, Type::Tensor);
			_value |= set(attr->second, "Tensor6", type, Type::Tensor6);
			_value |= set(attr->second, "Matrix", type, Type::Matrix);
			_value |= set(attr->second, "GlobalID", type, Type::GlobalID);
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "Center")) {
			_name = true;
			_value |= set(attr->second, "Node", center, Center::Node);
			_value |= set(attr->second, "Cell", center, Center::Cell);
			_value |= set(attr->second, "Grid", center, Center::Grid);
			_value |= set(attr->second, "Face", center, Center::Face);
			_value |= set(attr->second, "Edge", center, Center::Edge);
			_value |= set(attr->second, "Other", center, Center::Other);
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "ItemType")) {
			_name = true;
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "ElementFamily")) {
			_name = true;
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "ElementDegree")) {
			_name = true;
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "ElementCell")) {
			_name = true;
		}
		if (_name == false) {
			eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFAttribute.\n", attr->first.c_str(), attr->second.c_str());
		}
		if (_name == true && _value == false) {
			eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFAttribute.\n", attr->first.c_str());
		}
	}
	XDMFElement::parse(e);
}
