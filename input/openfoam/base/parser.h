#ifndef PARSER_H
#define PARSER_H

#include "tokenizer.h"
#include "../openfoam.h"
#include "../../loader.h"
#include <vector>

typedef std::vector< mesh::Point > Points;


template<typename T> ParseError* parse(Tokenizer &ts, T &value)
{
    return value.parse(ts);
}

template<> inline ParseError* parse(Tokenizer &ts, int &value)
{
    return ts.readInt(value);
}

template<> inline ParseError* parse(Tokenizer &ts, long &value)
{
    return ts.readLong(value);
}

template<> inline ParseError* parse(Tokenizer &ts, unsigned int &value)
{
    return ts.readInt(value);
}

template<> inline ParseError* parse(Tokenizer &ts, double &value)
{
    return ts.readDouble(value);
}

template<> inline ParseError* parse(Tokenizer &ts, bool &value)
{
    std::string s;
    PARSE_GUARD(ts.getIdentifier(s));

    if (s == "true")
    {
        value = true;
        return NULL;
    }

    if (s == "false")
    {
        value = false;
        return NULL;
    }

    return ts.makeError("'true' or 'false' expected");
}

template<> inline ParseError* parse(Tokenizer &ts, float &value)
{
    double d;
    PARSE_GUARD(ts.readDouble(d));
    value = d;
    return NULL;
}

template<> inline ParseError* parse(Tokenizer &ts, std::string &value)
{
    return ts.readIdentifier(value);
}

template<> inline ParseError* parse(Tokenizer &ts, mesh::Point &value)
{
    double v;
    PARSE_GUARD(ts.consumeChar('('));
    PARSE_GUARD(parse(ts, v));
    value.x=v;
    PARSE_GUARD(parse(ts, v));
    value.y=v;
    PARSE_GUARD(parse(ts, v));
    value.z=v;
    return ts.consumeChar(')');
}

template<typename T> ParseError* parse(Tokenizer &ts, T* &value)
{
    value = new T();
    return parse(ts, *value);
}

template<typename T> ParseError* parse(Tokenizer &ts, std::vector<T> &value)
{

    if (ts.isTokenIdentifier())
    {
        PARSE_GUARD(ts.nextToken());
    }
    if (ts.isTokenInt())
    {
        int i;
        ts.readInt(i);
        value.reserve(i);
    }
    PARSE_GUARD(ts.consumeChar('('));
    while (!ts.isTokenChar(')') && !ts.isTokenEnd())
    {
        T item;
        PARSE_GUARD(parse(ts, item));
        value.push_back(item);
    }
    PARSE_GUARD(ts.consumeChar(')'));
    return NULL;
}

#endif // PARSER_H
