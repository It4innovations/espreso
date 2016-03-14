#include "cellzone.h"

using namespace espreso::input;

CellZone::CellZone()
{

}

ParseError *CellZone::loadFromDictionary(Dictionary &dictionary)
{
    name = dictionary.getName();
    PARSE_GUARD(dictionary.readEntry("cellLabels", _elementIndexes));
    return NULL;
}

ParseError* espreso::input::parse(Tokenizer &ts, CellZone &cellZone)
{
    Dictionary dictionary;

    PARSE_GUARD(parse(ts, dictionary));
    PARSE_GUARD(cellZone.loadFromDictionary(dictionary));
    return NULL;
}


