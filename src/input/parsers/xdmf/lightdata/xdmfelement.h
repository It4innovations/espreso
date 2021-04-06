
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFELEMENT_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFELEMENT_H_

#include "basis/utilities/xml.h"
#include "basis/utilities/parser.h"

#include <vector>
#include <functional>

namespace espreso {

class XDMFAttribute;
class XDMFDataItem;
class XDMFDomain;
class XDMFGeometry;
class XDMFGrid;
class XDMFInformation;
class XDMFTime;
class XDMFTopology;

class XDMFElement {
public:
	enum class EType {
		Attribute,
		DataItem,
		Domain,
		Geometry,
		Grid,
		Information,
		Time,
		Topology
	};

	virtual void parse(XML::Element *e);
	virtual ~XDMFElement();

	XDMFElement* get(const std::string &name);
	void recurse(std::function<void(XDMFElement *e, EType type)> fnc);

	std::vector<XDMFAttribute*> attribute;
	std::vector<XDMFDataItem*> dataitem;
	std::vector<XDMFDomain*> domain;
	std::vector<XDMFGeometry*> geomery;
	std::vector<XDMFGrid*> grid;
	std::vector<XDMFInformation*> information;
	std::vector<XDMFTime*> time;
	std::vector<XDMFTopology*> topology;

protected:
	template <typename TType>
	bool set(const std::string &setting, const std::string &value, TType &parameter, TType pvalue)
	{
		if (StringCompare::caseInsensitiveEq(setting, value)) {
			parameter = pvalue;
			return true;
		}
		return false;
	}
};

}



#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFELEMENT_H_ */
