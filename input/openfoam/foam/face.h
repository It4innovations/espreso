#ifndef FACE_H
#define FACE_H

#include "../base/parser.h"
#include "../base/tokenizer.h"

class Face
{
public:
    Face();
    virtual ~Face();

    ParseError* parse(Tokenizer &ts);
    bool containsLine(int x, int y);
    ParseError* nextPoint(int x, int y, int &next);

    friend inline std::ostream& operator<<(std::ostream& os, const Face& obj)
    {
        // write obj to stream
        os<<obj.numberOfPoints<<"("<<obj.p[0]<<","<<obj.p[1]<<","<<obj.p[2];
        if (obj.numberOfPoints==4) {
        os<<","<<obj.p[3];
        }
        os<<")";
        return os;
    }

    int p[4];
    int numberOfPoints;

protected:
private:
};

typedef std::vector<Face> Faces;


#endif // FACE_H
