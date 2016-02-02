#ifndef FILESTREAM_H
#define FILESTREAM_H

#include <string>
#include <zlib.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "simplestream.h"

class FileStream : public SimpleStream
{
public:

    FileStream(const std::string &filename);
    ~FileStream();

    char read() {
        if (pointer < size) {
            return buffer[pointer++];
        } else {
            if (file.eof()) {
                return 0;
            }
            file.read(buffer, BUFFER_SIZE);
            size = file.gcount();
            pointer = 1;
            return buffer[0];
        }
    }

    ByteArray read(size_t size) {
        size_t r = this->size - pointer;
        if (r >= size) {
            size_t p = pointer;
            pointer += size;
            return ByteArray(&buffer[p], size);
        } else {
            char* tmp = new char[size];
            memcpy(tmp, &buffer[pointer], r * sizeof(char));
            size = 0;
            pointer = 0;
            file.read(&tmp[r],size - r);
            ByteArray ba(tmp, size);
            delete[] tmp;
            return ba;
        }
    }

    size_t tell() {
        return (size_t)file.tellg() - size + pointer;
    }

    void seek(size_t pos) {
        size_t p = file.tellg();
        if (pos < p && pos >= p - size) {
            pointer = pos - (p - size);
        } else {
            file.seekg(pos);
            file.read(buffer, BUFFER_SIZE);
            size = file.gcount();
            pointer = 0;
        }
    }

private:

    size_t pointer;
    size_t size;
    static const size_t BUFFER_SIZE = 8 * 1024 * 1024;
    char *buffer;
    std::ifstream file;
};

class GzFileStream : public SimpleStream
{
public:
    GzFileStream(const std::string &filename);

    ~GzFileStream() {
        gzclose(file);
    }

    char read() {
        int c = gzgetc(file);
        if (c == -1) {
            return 0;
        } else {
            return c;
        }
    }

    size_t tell() {
        return gztell(file);
    }

    void seek(size_t pos) {
        gzseek(file, pos, SEEK_SET);
    }

    ByteArray read(size_t size) {
        ByteArray ba;
        //ba.resize(size);
        gzread(file, ba.data(), size);
        return ba;
    }


protected:
    gzFile file;
};


#endif // FILESTREAM_H
