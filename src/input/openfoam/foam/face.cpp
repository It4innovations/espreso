#include "face.h"

using namespace espreso::input;

Face::Face()
{
    numberOfPoints =0;
    _index.first = NULL;
    _index.second = 0;
}

Face::~Face()
{
    //dtor
}

ParseError* Face::parse(Tokenizer &ts)
{
    PARSE_GUARD(ts.readeslocal(numberOfPoints));
    PARSE_GUARD(ts.consumeChar('('));
    PARSE_GUARD(ts.readeslocal(p[0]));
    PARSE_GUARD(ts.readeslocal(p[1]));
    PARSE_GUARD(ts.readeslocal(p[2]));
    if (numberOfPoints == 4)
    {
        PARSE_GUARD(ts.readeslocal(p[3]));
    }
    else if (numberOfPoints > 4)
    {
        return new ParseError("Face with more than 4 corners encountered.", "object: Face");
    }
    PARSE_GUARD(ts.consumeChar(')'));
    return NULL;
}

bool Face::containsLine(eslocal x, eslocal y)
{
    for(int i=0; i<numberOfPoints; i++)
    {
        if (p[i]==x)
        {
            if (p[(i+1)%numberOfPoints]==y) return true;
            if (p[(i-1+numberOfPoints)%numberOfPoints]==y) return true;
        }
    }
    return false;
}

ParseError* Face::nextPoint(eslocal x, eslocal y, eslocal &next)
{
    for(int i=0; i<numberOfPoints; i++)
    {
        if (p[i]==x)
        {
            if (p[(i+1)%numberOfPoints]==y)
            {
                next = p[(i+2)%numberOfPoints];
                return NULL;
            }
            if (p[(i-1+numberOfPoints)%numberOfPoints]==y)
            {
                next = p[(i-2+numberOfPoints)%numberOfPoints];
                return NULL;
            }
        }
    }
    std::stringstream ss;
    ss<<"No next point for line ("<<x<<","<<y<<") in Face: "<<this;
    return new ParseError(ss.str(), "Faces");
}

