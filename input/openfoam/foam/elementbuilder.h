#ifndef ELEMENTBUILDER_H
#define ELEMENTBUILDER_H

#include "face.h"

class ElementBuilder
{
public:
    ElementBuilder(Faces *faces);
    virtual ~ElementBuilder();

    void add(int index)
    {
        selectedFaces[numberOfFaces]=&(faces->at(index));
        numberOfFaces++;
    }

    friend inline std::ostream& operator<<(std::ostream& os, const ElementBuilder& obj)
    {
        // write obj to stream
        os<<obj.numberOfFaces<<"(";
        for(int i=0;i<obj.numberOfFaces;i++) {
            if (i!=0) os<<",";
            os<<*(obj.selectedFaces[i]);
        }
        os<<")";
        return os;
    }
    int getNumberOfFaces() { return numberOfFaces;}

    ParseError* createElement();


protected:
private:

    ParseError* nextPoint(Face *origin, eslocal x, eslocal y, eslocal &nextPoint);

    Face** selectedFaces;
    int numberOfFaces;
    Faces *faces;
};

#endif // ELEMENTBUILDER_H
