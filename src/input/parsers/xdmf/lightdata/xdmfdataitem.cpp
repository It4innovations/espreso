
#include "xdmfdataitem.h"
#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"

#include <sstream>

using namespace espreso;

XDMFDataItem::XDMFDataItem()
: itemtype(ItemType::Uniform), numbertype(NumberType::Float), precision(4), format(Format::XML), endian(Endian::Native), compression(Compression::Raw), seek(0)
{

}

void XDMFDataItem::parse(XML::Element *e)
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
        if (StringCompare::caseInsensitiveEq(attr->first, "Type") || StringCompare::caseInsensitiveEq(attr->first, "ItemType")) {
            _name = true;
            _value |= set(attr->second, "Uniform", itemtype, ItemType::Uniform);
            _value |= set(attr->second, "Collection", itemtype, ItemType::Collection);
            _value |= set(attr->second, "Tree", itemtype, ItemType::Tree);
            _value |= set(attr->second, "HyperSlab", itemtype, ItemType::HyperSlab);
            _value |= set(attr->second, "Coordinates", itemtype, ItemType::Coordinates);
            _value |= set(attr->second, "Function", itemtype, ItemType::Function);
        }
        if (StringCompare::caseInsensitiveEq(attr->first, "Dimensions")) {
            _name = _value = true;
            int dimension;
            std::stringstream ss(attr->second);
            while (ss.good()) {
                ss >> dimension;
                dimensions.push_back(dimension);
            }
        }
        if (StringCompare::caseInsensitiveEq(attr->first, "NumberType") || StringCompare::caseInsensitiveEq(attr->first, "DataType")) {
            _name = true;
            _value |= set(attr->second, "Float", numbertype, NumberType::Float);
            _value |= set(attr->second, "Int", numbertype, NumberType::Int);
            _value |= set(attr->second, "UInt", numbertype, NumberType::UInt);
            _value |= set(attr->second, "Char", numbertype, NumberType::Char);
            _value |= set(attr->second, "UChar", numbertype, NumberType::UChar);
        }
        if (StringCompare::caseInsensitiveEq(attr->first, "Precision")) {
            _name = _value = true;
            std::stringstream ss(attr->second);
            ss >> precision;
        }
        if (StringCompare::caseInsensitiveEq(attr->first, "Format")) {
            _name = true;
            _value |= set(attr->second, "XML", format, Format::XML);
            _value |= set(attr->second, "HDF", format, Format::HDF);
            _value |= set(attr->second, "Binary", format, Format::Binary);
        }
        if (StringCompare::caseInsensitiveEq(attr->first, "Endian")) {
            _name = true;
            _value |= set(attr->second, "Native", endian, Endian::Native);
            _value |= set(attr->second, "Big", endian, Endian::Big);
            _value |= set(attr->second, "Little", endian, Endian::Little);
        }
        if (StringCompare::caseInsensitiveEq(attr->first, "Compression")) {
            _name = true;
            _value |= set(attr->second, "Raw", compression, Compression::Raw);
            _value |= set(attr->second, "Zlib", compression, Compression::Zlib);
            _value |= set(attr->second, "BZip2", compression, Compression::BZip2);
        }
        if (_name == false) {
            eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFDataItem.\n", attr->first.c_str(), attr->second.c_str());
        }
        if (_name == true && _value == false) {
            eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFDataItem.\n", attr->first.c_str());
        }
    }
    data = e->value;
    XDMFElement::parse(e);
}
