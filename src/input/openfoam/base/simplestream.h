#ifndef SIMPLESTREAM_H
#define SIMPLESTREAM_H

#include <string>
#include "bytearray.h"

namespace espreso {
namespace input {

class SimpleStream
{
public:
    SimpleStream(const std::string &name) : name(name) {}
    virtual ~SimpleStream() {}

    const std::string &getName() {
        return name;
    }

    virtual char read() = 0;

protected:
    std::string name;
};


class ByteArrayStream : public SimpleStream
{
public:
    ByteArrayStream(const ByteArray &data, const std::string &name);

    char read()
    {
        if (position >= data.size()) {
            return 0;
        } else {
            return data.at(position++);
        }
    }

    ByteArray read(size_t size) {
        ByteArray ba(data.data() + position, size);
        position += size;
        return ba;
    }

    size_t tell() {
        return position;
    }

    void seek(size_t pos) {
        position = pos;
    }

protected:
    int position;
    ByteArray data;
};

}
}

#endif // SIMPLESTREAM_H
