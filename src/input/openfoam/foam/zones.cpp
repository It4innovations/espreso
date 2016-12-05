#include "zones.h"

using namespace espreso::input;

Zone::Zone(std::string containerIndexesName)
{
	this->containerIndexesName = containerIndexesName;
}

ParseError *Zone::loadFromDictionary(Dictionary &dictionary)
{
    name = dictionary.getName();
    PARSE_GUARD(dictionary.readEntry(containerIndexesName, _elementIndexes));
    return NULL;
}

ParseError* espreso::input::parse(Tokenizer &ts, Zone &zone)
{
    Dictionary dictionary;

    PARSE_GUARD(parse(ts, dictionary));
    PARSE_GUARD(zone.loadFromDictionary(dictionary));
    return NULL;
}

ParseError* espreso::input::parse(Tokenizer &ts, CellZone &cellZone)
{
    return parse(ts, (Zone&)cellZone);
}

ParseError* espreso::input::parse(Tokenizer &ts, FaceZone &faceZone)
{
    return parse(ts, (Zone&)faceZone);
}

ParseError* espreso::input::parse(Tokenizer &ts, PointZone &pointZone)
{
    return parse(ts, (Zone&)pointZone);
}
