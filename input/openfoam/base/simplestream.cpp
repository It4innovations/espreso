#include "simplestream.h"

ByteArrayStream::ByteArrayStream(const ByteArray &data, const std::string &name)
    : SimpleStream(name), position(0), data(data)
{

}

