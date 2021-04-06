
#include "xdmfgeometry.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

XDMFGeometry::XDMFGeometry()
: type(Type::XYZ)
{

}

void XDMFGeometry::parse(XML::Element *e)
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
		if (StringCompare::caseInsensitiveEq(attr->first, "Type") || StringCompare::caseInsensitiveEq(attr->first, "GeometryType")) {
			_name = true;
			_value |= set(attr->second, "XYZ", type, Type::XYZ);
			_value |= set(attr->second, "XY", type, Type::XY);
			_value |= set(attr->second, "X_Y_Z", type, Type::X_Y_Z);
			_value |= set(attr->second, "VxVyVzm", type, Type::VxVyVzm);
			_value |= set(attr->second, "Origin_DxDyDz", type, Type::Origin_DxDyDz);
			_value |= set(attr->second, "Origin_DxDy", type, Type::Origin_DxDy);
		}
		if (_name == false) {
			eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFGeometry.\n", attr->first.c_str(), attr->second.c_str());
		}
		if (_name == true && _value == false) {
			eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFGeometry.\n", attr->first.c_str());
		}
	}
	XDMFElement::parse(e);
}
