
#include "lightdata.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"

#include <sstream>

using namespace espreso;

LightData::LightData(const std::string &filename)
: version(3), filenane(filename), dir(utils::getFileDirectory(filename))
{
    XML xml;
    xml.load(filename);

    for (auto attr = xml.root.attributes.begin(); attr != xml.root.attributes.end(); ++attr) {
        bool _name = false, _value = false;
        if (StringCompare::caseInsensitiveEq(attr->first, "xmlns:xi")) {
            _name = _value = true;
        }
        if (StringCompare::caseInsensitiveEq(attr->first, "version")) {
            _name = _value = true;
            std::stringstream ss(attr->second);
            ss >> version;
        }
        if (_name == false) {
            eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFRoot.\n", attr->first.c_str(), attr->second.c_str());
        }
        if (_name == true && _value == false) {
            eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFRoot.\n", attr->first.c_str());
        }
    }
    XDMFElement::parse(&xml.root);
}
