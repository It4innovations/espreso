#ifndef ELEMENTBUILDER_H
#define ELEMENTBUILDER_H

#include "face.h"
#include <set>
#include <list>


class ElementBuilder
{
public:
    ElementBuilder();
    virtual ~ElementBuilder();

    void add(Face *face, bool owner)
    {
        //selectedFaces[numberOfFaces]=&(faces->at(index));
    	selectedFaces.push_back(std::pair<Face*, bool>(face, owner));
    }

    friend inline std::ostream& operator<<(std::ostream& os, const ElementBuilder& obj)
    {
        // write obj to stream
        os<<obj.selectedFaces.size()<<"(";
        bool first = true;
        for(std::list< std::pair<Face*, bool> >::const_iterator it = obj.selectedFaces.begin();
        				it != obj.selectedFaces.end(); ++it) {
            if (first) {
            	first=false;
            }else {
            	os<<",";
            }
            os<<*((*it).first)<<"-"<<(*it).second;
        }
        os<<")";
        return os;
    }
    size_t getNumberOfFaces() { return selectedFaces.size();}

    ParseError* createElement(std::vector<mesh::Element*> &elements);


protected:
private:

    ParseError* nextPoint(eslocal x, eslocal y, eslocal &nextPoint);

    std::list< std::pair<Face*, bool> > selectedFaces;
};

#endif // ELEMENTBUILDER_H
