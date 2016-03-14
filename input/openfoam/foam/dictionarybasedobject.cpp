#include "dictionarybasedobject.h"

using namespace espreso::input;

DictionaryBasedObject::DictionaryBasedObject()
{
}

DictionaryBasedObject::~DictionaryBasedObject()
{

}

ParseError *DictionaryBasedObject::parse(Tokenizer &ts)
{
    Dictionary dictionary;
    PARSE_GUARD(::parse(ts, dictionary));
    PARSE_GUARD(readDictionary(dictionary));
    return NULL;
}

void DictionaryBasedObject::write(TextStream &ts) const
{
    Dictionary dictionary;
    writeDictionary(dictionary);
    ::write(ts, dictionary);
}
