
#include "triangle6.h"

using namespace espreso;

Triangle6Generator::Triangle6Generator()
{
    subelements = 2;
    enodes = 6;
    code = Element::CODE::TRIANGLE6;
}

void Triangle6Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
    elements.push_back(indices[0]);
    elements.push_back(indices[2]);
    elements.push_back(indices[8]);
    elements.push_back(indices[1]);
    elements.push_back(indices[5]);
    elements.push_back(indices[4]);

    elements.push_back(indices[0]);
    elements.push_back(indices[8]);
    elements.push_back(indices[6]);
    elements.push_back(indices[4]);
    elements.push_back(indices[7]);
    elements.push_back(indices[3]);
}



