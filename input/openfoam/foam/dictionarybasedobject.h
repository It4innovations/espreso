#ifndef DIRECTORYOBJECT_H
#define DIRECTORYOBJECT_H

#include "dictionary.h"

class DictionaryBasedObject
{
public:
    DictionaryBasedObject();
    virtual ~DictionaryBasedObject();

    ParseError* parse(Tokenizer &ts);
    void write(TextStream &ts) const;    

    Dictionary* getDictionary() const {
        Dictionary *dictionary = new Dictionary();
        writeDictionary(*dictionary);
        return dictionary;
    }

    virtual ParseError* readDictionary(const Dictionary &dictionary) = 0;
    virtual void writeDictionary(Dictionary &dictionary) const = 0;
};

#endif // DIRECTORYOBJECT_H
