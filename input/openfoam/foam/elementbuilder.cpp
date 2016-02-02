#include "elementbuilder.h"
#include <set>

ElementBuilder::ElementBuilder(Faces *faces) : faces(faces)
{
    selectedFaces = new Face*[6];
    numberOfFaces=0;
}

ElementBuilder::~ElementBuilder()
{
    delete selectedFaces;
}

ParseError* ElementBuilder::createElement()
{

    std::set< int > coordinates;
    int numberOfSquares =0;
    int indicies[8];

    for(int i=0; i<numberOfFaces; i++)
    {
        Face *face=selectedFaces[i];
        coordinates.insert(face->p[0]);
        coordinates.insert(face->p[1]);
        coordinates.insert(face->p[2]);
        if (face->numberOfPoints==4)
        {
            coordinates.insert(face->p[3]);
            numberOfSquares++;
        }
    }
    if (coordinates.size()==4)
    {
        if (numberOfFaces!=4)
        {
            std::stringstream ss;
            ss<<"Element with 4 unique coordinates can not have "<<numberOfFaces<<" faces.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        if (numberOfSquares>0)
        {
            std::stringstream ss;
            ss<<"Element with 4 unique coordinates can not have a face with 4 points.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        Face *face = selectedFaces[0];
        coordinates.erase(face->p[0]);
        coordinates.erase(face->p[1]);
        coordinates.erase(face->p[2]);

        indicies[0] = face->p[0];
        indicies[1] = face->p[1];
        indicies[2] = face->p[2];
        indicies[3] = indicies[2];
        indicies[4] = *(coordinates.begin());
        //indicies[5] = indicies[4];
        //indicies[6] = indicies[4];
        //indicies[7] = indicies[4];
        //Tetrahedron4
        std::cout<<"Tetrahedron:";
        for(int i=0;i<5;i++) std::cout<<indicies[i]<<" ";
        std::cout<<"\n";
    }
    else if (coordinates.size()==5)
    {
        if (numberOfFaces!=5)
        {
            std::stringstream ss;
            ss<<"Element with 5 unique coordinates can not have "<<numberOfFaces<<" faces.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        if (numberOfSquares!=1)
        {
            std::stringstream ss;
            ss<<"Element with 5 unique coordinates must have 1 face with 4 points.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        Face *face;
        for(int i=0;i<numberOfFaces;i++) {
            Face *tmp = selectedFaces[i];
            if (tmp->numberOfPoints==4) {
                face = tmp;
                break;
            }
        }
        coordinates.erase(face->p[0]);
        coordinates.erase(face->p[1]);
        coordinates.erase(face->p[2]);
        coordinates.erase(face->p[3]);

        indicies[0] = face->p[0];
        indicies[1] = face->p[1];
        indicies[2] = face->p[2];
        indicies[3] = face->p[3];
        indicies[4] = *(coordinates.begin());

        //Pyramid5
        std::cout<<"Pyramid5:";
        for(int i=0;i<5;i++) std::cout<<indicies[i]<<" ";
        std::cout<<"\n";
    }
    else if (coordinates.size()==6)
    {
        if (numberOfFaces!=5)
        {
            std::stringstream ss;
            ss<<"Element with 6 unique coordinates can not have "<<numberOfFaces<<" faces.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        if (numberOfSquares!=3)
        {
            std::stringstream ss;
            ss<<"Element with 6 unique coordinates must have 3 faces with 4 points.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        Face *face;
        for(int i=0;i<numberOfFaces;i++) {
            Face *tmp = selectedFaces[i];
            if (tmp->numberOfPoints==3) {
                face = tmp;
                break;
            }
        }
        indicies[0] = face->p[0];
        indicies[1] = face->p[1];
        indicies[2] = face->p[2];
        indicies[3] = indicies[2];
        PARSE_GUARD(nextPoint(face, indicies[3], indicies[0],indicies[4]));
        PARSE_GUARD(nextPoint(face, indicies[0], indicies[1],indicies[5]));
        PARSE_GUARD(nextPoint(face, indicies[1], indicies[2],indicies[6]));
        indicies[7] = indicies[6];
        //Prism6
        std::cout<<"Prism:";
        for(int i=0;i<8;i++) std::cout<<indicies[i]<<" ";
        std::cout<<"\n";
    }
    else if (coordinates.size()==8)
    {
        if (numberOfFaces!=6)
        {
            std::stringstream ss;
            ss<<"Element with 8 unique coordinates can not have "<<numberOfFaces<<" faces.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        if (numberOfSquares!=6)
        {
            std::stringstream ss;
            ss<<"Element with 4 unique coordinates supports only faces with 4 points.";
            return new ParseError(ss.str(), "ElementBuilder");
        }
        indicies[0] = selectedFaces[0]->p[0];
        indicies[1] = selectedFaces[0]->p[1];
        indicies[2] = selectedFaces[0]->p[2];
        indicies[3] = selectedFaces[0]->p[3];

        PARSE_GUARD(nextPoint(selectedFaces[0], indicies[3], indicies[0],indicies[4]));
        PARSE_GUARD(nextPoint(selectedFaces[0], indicies[0], indicies[1],indicies[5]));
        PARSE_GUARD(nextPoint(selectedFaces[0], indicies[1], indicies[2],indicies[6]));
        PARSE_GUARD(nextPoint(selectedFaces[0], indicies[2], indicies[3],indicies[7]));
        //Hexahedron
        std::cout<<"Hexahedron: ";
        for(int i=0;i<8;i++) std::cout<<indicies[i]<<" ";
        std::cout<<"\n";
    }
    else
    {
        std::stringstream ss;
        ss<<"Element with "<<coordinates.size()<<" coordinates is not supported.";
        return new ParseError(ss.str(), "ElementBuilder");
    }
    return NULL;
}

ParseError* ElementBuilder::nextPoint(Face *origin, int x, int y, int &nextPoint) {
    for(int i=0;i<numberOfFaces;i++) {
        if (selectedFaces[i]!=origin) {
            if (selectedFaces[i]->containsLine(x,y)) {
                return selectedFaces[i]->nextPoint(x,y,nextPoint);
            }
        }
    }
    std::stringstream ss;
    ss<<"No next point for line ("<<x<<","<<y<<") in Element: "<<this;
    return new ParseError(ss.str(), "ElementBuilder");
}
