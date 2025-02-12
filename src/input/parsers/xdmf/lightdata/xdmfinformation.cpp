
#include "xdmfinformation.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

XDMFInformation::XDMFInformation()
{

}

void XDMFInformation::parse(XML::Element *e)
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
        if (StringCompare::caseInsensitiveEq(attr->first, "value")) {
            _name = _value = true;
            value = attr->second;
        }
        if (_name == false) {
            eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFInformation.\n", attr->first.c_str(), attr->second.c_str());
        }
        if (_name == true && _value == false) {
            eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFInformation.\n", attr->first.c_str());
        }
    }
    XDMFElement::parse(e);
}
