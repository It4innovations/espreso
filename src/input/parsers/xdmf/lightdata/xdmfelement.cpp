
#include <input/parsers/xdmf/lightdata/lightdata.h>
#include "xdmfelement.h"
#include "basis/utilities/parser.h"

#include "xdmfattribute.h"
#include "xdmfdataitem.h"
#include "xdmfdomain.h"
#include "xdmfgeometry.h"
#include "xdmfgrid.h"
#include "xdmfinformation.h"
#include "xdmftime.h"
#include "xdmftopology.h"

using namespace espreso;

void XDMFElement::parse(XML::Element *e)
{
    for (auto subelement = e->elements.begin(); subelement != e->elements.end(); ++subelement) {
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "Atrribute")) {
            attribute.push_back(new XDMFAttribute());
            attribute.back()->parse(*subelement);
        }
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "DataItem")) {
            dataitem.push_back(new XDMFDataItem());
            dataitem.back()->parse(*subelement);
        }
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "Domain")) {
            domain.push_back(new XDMFDomain());
            domain.back()->parse(*subelement);
        }
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "Geometry")) {
            geomery.push_back(new XDMFGeometry());
            geomery.back()->parse(*subelement);
        }
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "Grid")) {
            grid.push_back(new XDMFGrid());
            grid.back()->parse(*subelement);
        }
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "Information")) {
            information.push_back(new XDMFInformation());
            information.back()->parse(*subelement);
        }
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "Time")) {
            time.push_back(new XDMFTime());
            time.back()->parse(*subelement);
        }
        if (StringCompare::caseInsensitiveEq((*subelement)->name, "Topology")) {
            topology.push_back(new XDMFTopology());
            topology.back()->parse(*subelement);
        }
    }
}

XDMFElement* XDMFElement::get(const std::string &name)
{
    for (auto it = attribute.begin(); it != attribute.end(); ++it) {
        if (StringCompare::caseInsensitiveEq((*it)->name, name)) {
            return *it;
        }
    }
    for (auto it = dataitem.begin(); it != dataitem.end(); ++it) {
        if (StringCompare::caseInsensitiveEq((*it)->name, name)) {
            return *it;
        }
    }
    for (auto it = domain.begin(); it != domain.end(); ++it) {
        if (StringCompare::caseInsensitiveEq((*it)->name, name)) {
            return *it;
        }
    }
    for (auto it = geomery.begin(); it != geomery.end(); ++it) {
        if (StringCompare::caseInsensitiveEq((*it)->name, name)) {
            return *it;
        }
    }
    for (auto it = grid.begin(); it != grid.end(); ++it) {
        if (StringCompare::caseInsensitiveEq((*it)->name, name)) {
            return *it;
        }
    }
    for (auto it = information.begin(); it != information.end(); ++it) {
        if (StringCompare::caseInsensitiveEq((*it)->name, name)) {
            return *it;
        }
    }
    for (auto it = topology.begin(); it != topology.end(); ++it) {
        if (StringCompare::caseInsensitiveEq((*it)->name, name)) {
            return *it;
        }
    }
    return NULL;
}

void XDMFElement::recurse(std::function<void(XDMFElement *e, EType type)> fnc)
{
    for (auto it = attribute.begin(); it != attribute.end(); ++it) {
        fnc(*it, EType::Attribute); (*it)->recurse(fnc);
    }
    for (auto it = dataitem.begin(); it != dataitem.end(); ++it) {
        fnc(*it, EType::DataItem); (*it)->recurse(fnc);
    }
    for (auto it = domain.begin(); it != domain.end(); ++it) {
        fnc(*it, EType::Domain); (*it)->recurse(fnc);
    }
    for (auto it = geomery.begin(); it != geomery.end(); ++it) {
        fnc(*it, EType::Geometry); (*it)->recurse(fnc);
    }
    for (auto it = grid.begin(); it != grid.end(); ++it) {
        fnc(*it, EType::Grid); (*it)->recurse(fnc);
    }
    for (auto it = information.begin(); it != information.end(); ++it) {
        fnc(*it, EType::Information); (*it)->recurse(fnc);
    }
    for (auto it = time.begin(); it != time.end(); ++it) {
        fnc(*it, EType::Time); (*it)->recurse(fnc);
    }
    for (auto it = topology.begin(); it != topology.end(); ++it) {
        fnc(*it, EType::Topology); (*it)->recurse(fnc);
    }
}

XDMFElement::~XDMFElement()
{
    for (auto it = attribute.begin(); it != attribute.end(); ++it) { delete *it; }
    for (auto it = dataitem.begin(); it != dataitem.end(); ++it) { delete *it; }
    for (auto it = domain.begin(); it != domain.end(); ++it) { delete *it; }
    for (auto it = geomery.begin(); it != geomery.end(); ++it) { delete *it; }
    for (auto it = grid.begin(); it != grid.end(); ++it) { delete *it; }
    for (auto it = information.begin(); it != information.end(); ++it) { delete *it; }
    for (auto it = time.begin(); it != time.end(); ++it) { delete *it; }
    for (auto it = topology.begin(); it != topology.end(); ++it) { delete *it; }
}
