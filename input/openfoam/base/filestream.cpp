#include "filestream.h"

using namespace espreso::input;

FileStream::FileStream(const std::string &filename)
    : SimpleStream(filename), file(filename.c_str())
{
    if (!file.is_open()){
        fprintf(stderr,"FileStream - can not open %s: ",filename.c_str());
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
        fprintf(stderr,"GzFileStream open failed: %s",filename.c_str());
        exit(0);
    }
}
