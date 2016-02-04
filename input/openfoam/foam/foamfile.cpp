
#include "foamfile.h"
#include "../base/filestream.h"
#include "../base/tokenstream.h"

bool FoamFile::fileExists(const char *filename)
{
    std::ifstream ifile(filename);
    return ifile;
}

FoamFile::FoamFile(const std::string &filename)
{
    printf("Loading %s ...\n", filename.c_str());

    std::string gzFilename = filename + ".gz";
    if (fileExists(gzFilename.c_str())) {
        GzFileStream *s = new GzFileStream(gzFilename);
        tokenizer = new TokenStream<GzFileStream>(s);
        stream = s;
    } else {
        FileStream *s = new FileStream(filename);
        tokenizer = new TokenStream<FileStream>(s);
        stream = s;
    }
    parseHeader();
}

std::string FoamFile::checkFileType(const std::string &filename)
{
    std::string gzFilename = filename + ".gz";
    if (fileExists(gzFilename.c_str())) {
        return gzFilename;
    } else if (fileExists(filename.c_str())) {
        return filename;
    }
    return "";
}

FoamFile::~FoamFile()
{
    delete tokenizer;
}

void FoamFile::parseHeader()
{
    std::string ident;
    ParseError *err = tokenizer->readIdentifier(ident);
    if (err) {
        fprintf(stderr, "parseHeader: %s\n", err->getMessage().c_str());
        exit(1);
    }
    if (ident != "FoamFile") {
        fprintf(stderr, "Invalid foam file\n");
        exit(1);
    }
    do {
        tokenizer->nextToken();
        if (tokenizer->failIfEnd() != NULL) {
            fprintf(stderr, "Invalid FoamFile\n");
            exit(1);
        }
    } while (!tokenizer->isTokenChar('}'));
    tokenizer->nextToken();
}
