#include "filestream.h"

#include "../../../basis/logging/logging.h"

using namespace espreso::input;

FileStream::FileStream(const std::string &filename)
    : SimpleStream(filename), file(filename.c_str())
{
    if (!file.is_open()){
        ESINFO(ERROR) << "FileStream - can not open: " << filename;
        exit(0);
    }
    buffer = new char[BUFFER_SIZE];
    pointer = 0;
    size = 0;
}

FileStream::~FileStream()
{
    delete [] buffer;
}


GzFileStream::GzFileStream(const std::string &filename)
    : SimpleStream(filename)
{
    file = gzopen(filename.c_str(), "r");
    if (file == NULL) {
        ESINFO(ERROR) << "GzFileStream open failed: " << filename;
        exit(0);
    }
}
