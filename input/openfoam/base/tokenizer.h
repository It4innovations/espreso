#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <string>
#include <iostream>
#include <stdlib.h>
#include <sstream>

#include "parseerror.h"
#include "simplestream.h"

namespace espreso {
namespace input {

enum TokenType {
    TOKEN_INT = 0,
    TOKEN_DOUBLE,
    TOKEN_CHAR,
    TOKEN_IDENT,
    TOKEN_END
};

inline char is_white_space(char c)
{
    return (c == ' ') || (c == '\t') || (c == '\n');
}

class Tokenizer
{
public:

    virtual ParseError *nextToken() = 0;
    virtual ~Tokenizer() {}

    virtual size_t tell() = 0;
    virtual void readRawText(size_t from, size_t to, ByteArray &out) = 0;

    bool isTokenEnd() {
        return token == TOKEN_END;
    }

    bool isTokenInt() {
        return token == TOKEN_INT;
    }

    bool isTokenIdentifier() {
        return token == TOKEN_IDENT;
    }

    bool isTokenChar(char c) {
        return token == TOKEN_CHAR && intValue == c;
    }

    ParseError *failIfEnd()
    {
        if (isTokenEnd()) {
            std::cerr << "Unexpected end of file." << std::endl;
            exit(1);
        }
        return NULL;
    }

    ParseError* getIdentifier(std::string &value) {
        PARSE_GUARD(expect(TOKEN_IDENT));
        value = stringValue;
        return NULL;
    }

    ParseError* expect(TokenType t) {
        if (token != t) {
            std::stringstream ss;
            ss << "Expected token " << t << ", got" << debugString();
            return makeError(ss.str());
        }
        return NULL;
    }

    ParseError *readIdentifier(std::string &text)
    {
        PARSE_GUARD(expect(TOKEN_IDENT));
        text = stringValue;
        PARSE_GUARD(nextToken());
        return NULL;
    }

    ParseError *readInt(int &value)
    {
        PARSE_GUARD(expect(TOKEN_INT));
        value = intValue;
        PARSE_GUARD(nextToken());
        return NULL;
    }

    ParseError *readLong(long &value)
    {
        PARSE_GUARD(expect(TOKEN_INT));
        value = intValue;
        PARSE_GUARD(nextToken());
        return NULL;
    }


    ParseError *readInt(unsigned int &value)
    {
        PARSE_GUARD(expect(TOKEN_INT));
        value = intValue;
        PARSE_GUARD(nextToken());
        return NULL;
    }

    ParseError *readDouble(double &value)
    {
        if (token == TOKEN_INT) {
            value = intValue;
        } else if (token == TOKEN_DOUBLE) {
            value = doubleValue;
        } else {
            return makeError("Expected double");
        }
        PARSE_GUARD(nextToken());
        return NULL;
    }

    ParseError *consumeChar(char c)
    {
        if (!isTokenChar(c)) {
            std::stringstream ss;
            ss<<"Expected char "<< c << debugString();
            return makeError(ss.str());
        }
        nextToken();
        return NULL;
    }

    ParseError *consumeIdentifier(const std::string &value) {
        if (isTokenIdentifier() && stringValue == value) {
            nextToken();
            return NULL;
        }
        return makeError("Identifier '" + value + "' expected");
    }

    ParseError *skipInt() {
        if (isTokenInt()) {
            return nextToken();
        } else {
            return NULL;
        }
    }

    std::string debugString()
    {
        std::stringstream ss;
        switch (token) {
        case TOKEN_END:
            return "EOF";
        case TOKEN_CHAR:
            ss << "Char " << static_cast<char>(intValue);
            return ss.str();
        case TOKEN_INT:
            ss<<"Int " << intValue;
            return ss.str();
        case TOKEN_DOUBLE:
            ss<<"Double "<<doubleValue;
            return ss.str();
        case TOKEN_IDENT:
            ss << "Ident " << stringValue;
            return ss.str();
        default:
            ss<< "Invalid token " << token;
            return ss.str();
        }
    }

    ParseError *makeError(const std::string &message) {
        return new ParseError(message, getSource());
    }

    virtual std::string getSource() = 0;
protected:

    TokenType token;
    std::string stringValue;
    long int intValue;
    double doubleValue;
};

}
}

#endif // TOKENIZER_H
