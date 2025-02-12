
#ifndef SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFDATAITEM_H_
#define SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFDATAITEM_H_

#include "xdmfelement.h"
#include <string>
#include <vector>

namespace espreso {

class XDMFDataItem: public XDMFElement {
public:
    enum class ItemType    { Uniform, Collection, Tree, HyperSlab, Coordinates, Function };
    enum class NumberType  { Float, Int, UInt, Char, UChar };
    enum class Format      { XML, HDF, Binary };
    enum class Endian      { Native, Big, Little };
    enum class Compression { Raw, Zlib, BZip2 };

    std::string name;
    std::string reference;
    ItemType itemtype;
    std::vector<esint> dimensions;
    NumberType numbertype;
    int precision;
    Format format;
    Endian endian;
    Compression compression;
    esint seek;
    std::string data;

    XDMFDataItem();
    void parse(XML::Element *e);
};

}



#endif /* SRC_INPUT_FORMATS_XDMF_LIGHTDATA_XDMFDATAITEM_H_ */
