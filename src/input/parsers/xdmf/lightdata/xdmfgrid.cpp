
#include "xdmfgrid.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

XDMFGrid::XDMFGrid()
: type(Type::Uniform), collectiontype(CollectionType::Spatial), section(Section::DataItem)
{

}

void XDMFGrid::parse(XML::Element *e)
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
		if (StringCompare::caseInsensitiveEq(attr->first, "Type") || StringCompare::caseInsensitiveEq(attr->first, "GridType")) {
			_name = true;
			_value |= set(attr->second, "Uniform", type, Type::Uniform);
			_value |= set(attr->second, "Collection", type, Type::Collection);
			_value |= set(attr->second, "Tree", type, Type::Tree);
			_value |= set(attr->second, "Subset", type, Type::Subset);
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "CollectionType")) {
			_name = true;
			_value |= set(attr->second, "Spatial", collectiontype, CollectionType::Spatial);
			_value |= set(attr->second, "Temporal", collectiontype, CollectionType::Temporal);
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "Section")) {
			_name = true;
			_value |= set(attr->second, "DataItem", section, Section::DataItem);
			_value |= set(attr->second, "All", section, Section::All);
		}
		if (_name == false) {
			eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFGrid.\n", attr->first.c_str(), attr->second.c_str());
		}
		if (_name == true && _value == false) {
			eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFGrid.\n", attr->first.c_str());
		}
	}
	XDMFElement::parse(e);
}
