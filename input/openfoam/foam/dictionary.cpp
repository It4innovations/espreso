#include "dictionary.h"

Dictionary::~Dictionary()
{
    for (std::vector< Dictionary* >::iterator it = dictionaries.begin() ; it != dictionaries.end(); ++it) {
        delete *it;
    }
}

void Dictionary::write(TextStream &ts) const
{
    ts << name << "\n";
    ts.beginLine();
    ts << "{\n";
    ts.incIndent();
    for (size_t i = 0; i < entries.size(); i++) {
        ts.beginLine();
        ts << entries[i].first << " " << entries[i].second << ";\n";
    }
    for (size_t i = 0; i < dictionaries.size(); i++) {
        ts.beginLine();
        dictionaries[i]->write(ts);
        ts << "\n";
    }
    ts.decIndent();
    ts.beginLine();
    ts << "}\n";
    ts.beginLine();
}

const ByteArray *Dictionary::lookupConst(const std::string &entryName) const
{
    for (size_t i = 0; i < entries.size(); i++) {
        if (entries[i].first == entryName) {
            return &entries[i].second;
        }
    }
    return NULL;
}

ByteArray *Dictionary::lookup(const std::string &entryName)
{
    for (size_t i = 0; i < entries.size(); i++) {
        if (entries[i].first == entryName) {
            return &entries[i].second;
        }
    }
    return NULL;
}

ParseError *Dictionary::parseEntries(Tokenizer &ts)
{
    while (!ts.isTokenChar('}') && !ts.isTokenEnd()) {
        DictionaryEntry entry;
        PARSE_GUARD(ts.getIdentifier(entry.first));

        size_t start = ts.tell();
        int level = 0;
        ts.nextToken();

        if (ts.isTokenChar('{')) {
            // SubDir
            Dictionary *dictionary = new Dictionary(entry.first);
            dictionary->parse(ts);
            dictionaries.push_back(dictionary);
            continue;
        }

        // Entry
        for(;;) {
            if (ts.isTokenEnd()) {
                return ts.makeError("Unexpected end of file");
            }

            if (ts.isTokenChar(';') && level == 0) {
                break;
            }

            if (ts.isTokenChar('{')) {
                level += 1;
            }

            if (ts.isTokenChar('}')) {
                level -= 1;
                if (level < 0) {
                    return ts.makeError("Unmatched parenthesis");
                }
            }

            PARSE_GUARD(ts.nextToken());
        }
        size_t end = ts.tell();
        ts.readRawText(start, end, entry.second);
        PARSE_GUARD(ts.nextToken());
        entries.push_back(entry);
    }
    return NULL;
}

ParseError *Dictionary::parse(Tokenizer &ts)
{
    source = ts.getSource();
    if (ts.isTokenIdentifier()) {
        PARSE_GUARD(ts.readIdentifier(name));
    }
    PARSE_GUARD(ts.consumeChar('{'));
    PARSE_GUARD(parseEntries(ts));
    PARSE_GUARD(ts.consumeChar('}'));
    if (ts.isTokenChar(';')) {
        PARSE_GUARD(ts.nextToken());
    }
    return NULL;
}

ParseError *Dictionary::parseTopLevel(Tokenizer &ts)
{
    source = ts.getSource();
    return parseEntries(ts);
}

void Dictionary::setSubDict(Dictionary *dictionary)
{
    for (size_t i = 0; i < dictionaries.size(); i++) {
        if (dictionaries[i]->getName() == name) {
            delete dictionaries[i];
            dictionaries[i] = dictionary;
        }
    }
    dictionaries.push_back(dictionary);
}

Dictionary *Dictionary::getSubDict(const std::string &name)
{
    for (size_t i = 0; i < dictionaries.size(); i++) {
        if (dictionaries[i]->getName() == name) {
            return dictionaries[i];
        }
    }
    return NULL;
}

std::vector< std::string > Dictionary::getEntryNames()
{
    std::vector< std::string > names;
    for (size_t i = 0; i < entries.size(); i++) {
        names.push_back(entries[i].first);
    }
    return names;
}
