#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <vector>
#include <utility>
#include "../base/tokenstream.h"
#include "../base/simplestream.h"
#include "../base/parser.h"
#include "../base/textstream.h"

typedef std::pair<std::string, ByteArray> DictionaryEntry;


class Dictionary
{
public:
    Dictionary(const std::string &name = "") : name(name) {}
    ~Dictionary();

    const std::string & getName() const {
        return name;
    }

    void setName(const std::string &name) {
        this->name = name;
    }

    ParseError* parse(Tokenizer &ts);
    ParseError* parseTopLevel(Tokenizer &ts);

    TokenStream<ByteArrayStream>* getTokenStream(const std::string &entryName) {
        const ByteArray *text = lookup(entryName);
        if (text == NULL) {
            return NULL;
        }
        ByteArrayStream *stream = new ByteArrayStream(*text, source);
        return new TokenStream<ByteArrayStream>(stream);
    }


    template<typename T>
    ParseError* readEntry(const std::string &entryName, T &value) const {
        const ByteArray *text = lookupConst(entryName);
        if (text == NULL) {
            return new ParseError("Dictionary has no entry "+entryName, name);
        }
        return parse(text, entryName, value);
    }

    template<typename T>
    ParseError* readEntry(const std::string &entryName, T &value, const T &default_value) const {
        const ByteArray *text = lookupConst(entryName);
        if (text == NULL) {
            value = default_value;
            return NULL;
        }
        return parse(text, entryName, value);
    }

    /*template<typename T>
    void setEntry(const std::string &entryName, const T &value) {
        std::string s;
        TextStream ts(&s);
        ::write(ts, value);
        ByteArray *r = lookup(entryName);
        if (r != NULL) {
            *r = new ByteArray(s.c_str(),s.size());
        } else  {
            entries.append(DictionaryEntry(entryName, s.toLocal8Bit()));
        }
    }*/

    void setSubDict(Dictionary *dictionary);
    Dictionary* getSubDict(const std::string &name);

    template<typename T>
    ParseError *readAllEntries(std::vector<T> &values) const {
        values.clear();
        for (int i = 0; i < entries.size(); i++) {;
            values.push_back(T());
            PARSE_GUARD(parse(entries[i].second, entries[i].first, values[i]));
        }
        return NULL;
    }

    std::vector< std::string > getEntryNames();

    const std::vector<Dictionary*> & getSubDicts() {
        return dictionaries;
    }

    void write(TextStream &ts) const;
    const ByteArray *lookupConst(const std::string &entryName) const;
    ByteArray *lookup(const std::string &entryName);

protected:

    ParseError* parseEntries(Tokenizer &ts);

    template<typename T>
    ParseError* parse(const ByteArray *text, const std::string &entryName, T &value) const {
        std::stringstream ss;
        ss << source <<"["<< entryName << "]";
        std::string s = ss.str();
        if (text == NULL) {
            return new ParseError("Dictionary '" + entryName + "' entry not found", s);
        }
        ByteArrayStream *stream = new ByteArrayStream(*text, s);
        TokenStream<ByteArrayStream> ts(stream);
        return ::parse(ts, value);
    }

    std::vector<DictionaryEntry> entries;
    std::vector<Dictionary*> dictionaries;
    std::string name;
    std::string source;
};

#endif // DICTIONARY_H
