
#include "xdmftopology.h"
#include "esinfo/eslog.hpp"

#include <sstream>

using namespace espreso;

XDMFTopology::XDMFTopology()
: type(Type::Mixed), nodeperelement(0), numberofelement(0), dimension(0), order(0)
{

}

void XDMFTopology::parse(XML::Element *e)
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
		if (StringCompare::caseInsensitiveEq(attr->first, "Type") || StringCompare::caseInsensitiveEq(attr->first, "TopologyType")) {
			_name = true;
			_value |= set(attr->second, "Polyvertex", type, Type::Polyvertex);
			_value |= set(attr->second, "Polyline", type, Type::Polyline);
			_value |= set(attr->second, "Polygon", type, Type::Polygon);
			_value |= set(attr->second, "Triangle", type, Type::Triangle);
			_value |= set(attr->second, "Quadrilateral", type, Type::Quadrilateral);
			_value |= set(attr->second, "Tetrahedron", type, Type::Tetrahedron);
			_value |= set(attr->second, "Pyramid", type, Type::Pyramid);
			_value |= set(attr->second, "Wedge", type, Type::Wedge);
			_value |= set(attr->second, "Hexahedron", type, Type::Hexahedron);
			_value |= set(attr->second, "Quadrilateral", type, Type::Quadrilateral);
			_value |= set(attr->second, "Edge_3", type, Type::Edge_3);
			_value |= set(attr->second, "Tri_6", type, Type::Triangle_6);
			_value |= set(attr->second, "Quad_8", type, Type::Quadrilateral_8);
			_value |= set(attr->second, "Tet_10", type, Type::Tetrahedron_10);
			_value |= set(attr->second, "Pyramid_13", type, Type::Pyramid_13);
			_value |= set(attr->second, "Wedge_15", type, Type::Wedge_15);
			_value |= set(attr->second, "Hex_20", type, Type::Hexahedron_20);
			_value |= set(attr->second, "Mixed", type, Type::Mixed);
			_value |= set(attr->second, "2DSMesh", type, Type::SMesh2D);
			_value |= set(attr->second, "2DRectMesh", type, Type::RectMesh2D);
			_value |= set(attr->second, "2DCoRectMesh", type, Type::CoRectMesh2D);
			_value |= set(attr->second, "3DSMesh", type, Type::SMesh3D);
			_value |= set(attr->second, "3DRectMesh", type, Type::RectMesh3D);
			_value |= set(attr->second, "3DCoRectMesh", type, Type::CoRectMesh3D);
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "NodesPerElement")) {
			_name = _value = true;
			std::stringstream ss(attr->second);
			ss >> nodeperelement;
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "NumberOfElements") || StringCompare::caseInsensitiveEq(attr->first, "Dimensions")) {
			_name = _value = true;
			std::stringstream ss(attr->second);
			ss >> numberofelement;
		}
		if (StringCompare::caseInsensitiveEq(attr->first, "Order")) {
			_name = _value = true;
			std::stringstream ss(attr->second);
			ss >> order;
		}
		if (_name == false) {
			eslog::warning("XDMF Reader: unknown attribute '%s=%s' in class XDMFTopology.\n", attr->first.c_str(), attr->second.c_str());
		}
		if (_name == true && _value == false) {
			eslog::warning("XDMF Reader: unknown value of attribute '%s' in class XDMFTopology.\n", attr->first.c_str());
		}
	}
	XDMFElement::parse(e);
}
