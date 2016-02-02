#ifndef FOAMFILE_H
#define FOAMFILE_H

#include <string>
#include <stdlib.h>
#include "../base/tokenizer.h"
#include "../base/simplestream.h"

class FoamFile
{
public:
    FoamFile(const std::string &filename);
    ~FoamFile();

    Tokenizer& getTokenizer() {
        return *tokenizer;
    }

    static std::string checkFileType(const std::string &filename);

    static bool exists(const std::string &filename) {
        return checkFileType(filename) != "";
    }

protected:
    void parseHeader();
    static bool fileExists(const char *filename);

    SimpleStream *stream;
    Tokenizer *tokenizer;
};

#endif // FOAMFILE_H
