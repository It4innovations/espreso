#include "bytearray.h"
#include <string.h>

ByteArray::ByteArray(const char* data, int size)
{
    _size = size;
    _data = new char[_size];
    memcpy(_data, data, size * sizeof(char));
}

ByteArray::ByteArray(const ByteArray &b)
{
    _size = b.size();
    _data = new char[_size];
    memcpy(_data, b._data, _size * sizeof(char));
}


ByteArray::ByteArray()
{
    _data = NULL;
    _size = 0;
}

ByteArray::~ByteArray()
{
    delete[] _data;
}

ByteArray & ByteArray::operator= (const ByteArray &b)
    {
        if (this != &b)
        {
            _size = b.size();
            if (_data!=NULL) delete [] _data;
            _data = new char[_size];
            memcpy(_data, b._data, _size * sizeof(char));
        }
        return *this;
    }


