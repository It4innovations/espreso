#include "simplestream.h"

using namespace espreso::input;

ByteArrayStream::ByteArrayStream(const ByteArray &data, const std::string &name)
    : SimpleStream(name), position(0), data(data)
{

}

