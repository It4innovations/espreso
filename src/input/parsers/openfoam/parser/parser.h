
#ifndef SRC_INPUT_OPENFOAM_PARSER_PARSER_H_
#define SRC_INPUT_OPENFOAM_PARSER_PARSER_H_

#include <cstdlib>
#include <string>

namespace espreso {

struct OpenFOAMParser {

    const char *begin, *end;

    OpenFOAMParser(const char *begin, const char *end);

    double readDouble(const char* &c)
    {
        char *endptr;
        double value = strtod(c, &endptr);
        c = endptr;
        return value;
    }

    double readInteger(const char* &c)
    {
        char *endptr;
        esint value = strtol(c, &endptr, 10);
        c = endptr;
        return value;
    }

    bool isEmpty(const char* &c)
    {
        return
                *c == ' ' ||
                *c == '\n' ||
                *c == '\r' ||
                *c == '\t';
    }

    void toNext(const char* &c)
    {
        while (!isEmpty(c)) {
            ++c;
        }
    }

    std::string readString(const char* &c)
    {
        while (isEmpty(c)) { ++c; }
        const char *s = c;
        while (!isEmpty(c)) { ++c; }
        return std::string(s, c);
    }
};

struct OpenFOAMSeparateParser: public OpenFOAMParser {

    OpenFOAMSeparateParser(const char *begin, const char *end);
};

struct OpenFOAMCollectiveParser: public OpenFOAMParser {


    OpenFOAMCollectiveParser(const char *begin, const char *end);
};

}




#endif /* SRC_INPUT_OPENFOAM_PARSER_PARSER_H_ */
