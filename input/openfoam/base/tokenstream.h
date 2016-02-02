#ifndef TOKENSTREAM_H
#define TOKENSTREAM_H

#include "tokenizer.h"
#include <sstream>
#include <string>

template<typename TStream>
class TokenStream : public Tokenizer
{
public:

    TokenStream(TStream *stream) : stream(stream)
    {
        next();
        nextToken();
    }

    ~TokenStream() {
        delete stream;
    }

    ParseError *nextToken()
    {
        for (;;) {
            skipWhiteSpaces();

            if ((currentChar >= '0' && currentChar <= '9') || currentChar == '-') {
                return _readNumber();
            }

            switch (currentChar) {
            case '{': // pass through
            case '}': // pass through
            case '(': // pass through
            case ')': // pass through
            case ';':
                intValue = currentChar;
                next();
                token = TOKEN_CHAR;
                return NULL;
            case '/': {
                next();
                if (currentChar == '/') { // One line comment
                    find('\n');
                    continue;
                } else if (currentChar == '*') {
                    next();
                    do {
                        if (find('*') == 0) {
                            return makeError("Unclosed comment");
                        }
                        next();
                    } while(currentChar != '/');
                    next();
                    continue;
                } else {
                    return makeError("Invalid comment type");
                }
            }
            case 0:
                token = TOKEN_END;
                return NULL;
            }

            return _readIdentifier();
        }
        return NULL;
    }

    virtual size_t tell() {
        return stream->tell();
    }

    virtual void readRawText(size_t from, size_t to, ByteArray &out) {
        size_t current = stream->tell();
        stream->seek(from);
        out = stream->read(to - from);
        stream->seek(current);
    }

    std::string getSource() {
        return stream->getName();
    }

protected:

    void next() {
        currentChar = stream->read();
    }

    void skipWhiteSpaces()
    {
        while (is_white_space(currentChar)) {
            next();
        }
    }

    char find(char c)
    {
        while(c != currentChar && currentChar != 0) {
            next();
        }
        return currentChar;
    }

    ParseError* _readIdentifier()
    {
        char buffer[201];
        int i = 0;
        do {
            buffer[i++] = currentChar;
            if (i > 200) {
                return makeError("Too long identifier");
            }
            next();
        } while(currentChar != 0 &&
                currentChar != ';' &&
                currentChar != ')' &&
                currentChar != '}' &&
                !is_white_space(currentChar));
        buffer[i] = 0;
        stringValue = std::string(buffer);
        token = TOKEN_IDENT;
        return NULL;
    }

    ParseError* _readNumber()
    {
        char buffer[100];
        bool integer_value = true;
        int i;
        for(i = 0; i < 100; i++) {
            if ((currentChar < '0' || currentChar > '9') && currentChar != '-' && currentChar != '+') {
                if (currentChar != '.' && currentChar != 'e' && currentChar != 'E') {
                    break;
                }
                integer_value = false;
            }
            buffer[i] = currentChar;
            next();
        }

        if (i == 100) {
            return makeError("Invalid numeric value");
        }

        buffer[i] = 0;

        char *end;
        if (!integer_value) {
            doubleValue = strtod(buffer, &end);
            if (*end != 0) {
                return makeError("Invalid floating point value");
            }
            token = TOKEN_DOUBLE;
            return NULL;
        } else {
            intValue = strtol(buffer, &end, 10);
            if (*end != 0) {
                return makeError("Invalid integer value");
            }
            token = TOKEN_INT;
            return NULL;
        }
    }

    char currentChar;
    TStream *stream;

};


#endif // TOKENSTREAM_H
