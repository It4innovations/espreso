#ifndef FACE_H
#define FACE_H

#include "../base/parser.h"
#include "../base/tokenizer.h"
#include "../../loader.h"

class Face
{
public:
    Face();
    virtual ~Face();

    ParseError* parse(Tokenizer &ts);
    bool containsLine(eslocal x, eslocal y);
    ParseError* nextPoint(eslocal x, eslocal y, eslocal &next);
    void setFaceIndex(mesh::Element *element, unsigned char index) {
    	_index = std::pair<mesh::Element*, unsigned char>(element, index);
    }
    mesh::FaceIndex &getFaceIndex() {
    	return _index;
    }

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

    eslocal p[4];
    int numberOfPoints;

protected:
private:
    mesh::FaceIndex _index;

};

#endif // FACE_H
