#ifndef PARSEERROR_H
#define PARSEERROR_H

#include <string>
#include <stdio.h>
#include "esbasis.h"

#define PARSE_GUARD(x) { ParseError *_parse_err = (x); if (_parse_err != NULL) { return _parse_err; } }

namespace espreso {
namespace input {

class ParseError
{
public:
    ParseError(const std::string &message, const std::string &source) {
        std::stringstream ss;
        ss<<source<<": "<<message;
        this->message = ss.str();
    }

    const std::string& getMessage() { return message; }

    void print() {
        ESINFO(ERROR) << "ParseError: " << message;
    }

protected:
    std::string message;
};

inline void catchParseError(ParseError *parseError) {
    if (parseError) {
        parseError->print();
        exit(1);
    }
}

}
}


#endif // PARSEERROR_H
