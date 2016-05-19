#ifndef BYTEARRAY_H
#define BYTEARRAY_H

#include <ostream>

namespace espreso {
namespace input {

class ByteArray
{
public:
    ByteArray(const char* data, int size);
    ByteArray(const ByteArray &b);
    ByteArray();
    virtual ~ByteArray();

    ByteArray & operator= (const ByteArray &b);

    inline   int size() const
    {
        return _size;
    }
    inline char * data()
    {
        return _data;
    }
    inline char at(int index) const
    {
        return _data[index];
    }

friend inline std::ostream& operator<<(std::ostream& os, const ByteArray& obj)
{
        // write obj to stream
        os.write(obj._data, obj._size);
        return os;
}


protected:
private:
    int _size;
    char* _data;
};

}
}


#endif // BYTEARRAY_H
