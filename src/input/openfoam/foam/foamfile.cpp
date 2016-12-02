
#include "foamfile.h"
#include "../base/filestream.h"
#include "../base/tokenstream.h"

using namespace espreso::input;

bool FoamFile::fileExists(const char *filename)
{
    std::ifstream ifile(filename);
    return ifile.good();
}

FoamFile::FoamFile(const std::string &filename)
{
    std::string gzFilename = filename + ".gz";

    if (fileExists(gzFilename.c_str())) {
        GzFileStream *s = new GzFileStream(gzFilename);
        tokenizer = new TokenStream<GzFileStream>(s);
        stream = s;
        ESINFO(OVERVIEW) << "Reading OpenFoam file: "<<gzFilename;

    } else {
        FileStream *s = new FileStream(filename);
        tokenizer = new TokenStream<FileStream>(s);
        stream = s;
        ESINFO(OVERVIEW) << "Reading OpenFoam file: "<<filename;
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
        ESINFO(ERROR) << "parseHeader: " << err->getMessage();
        exit(1);
    }
    if (ident != "FoamFile") {
        ESINFO(ERROR) << "Invalid foam file";
        exit(1);
    }
    do {
        tokenizer->nextToken();
        if (tokenizer->failIfEnd() != NULL) {
            ESINFO(ERROR) << "Invalid FoamFile";
            exit(1);
        }
    } while (!tokenizer->isTokenChar('}'));
    tokenizer->nextToken();
}
