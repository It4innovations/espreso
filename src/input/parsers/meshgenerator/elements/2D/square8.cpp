
#include "square8.h"

using namespace espreso;

Square8Generator::Square8Generator()
{
    subelements = 1;
    enodes = 8;
    code = Element::CODE::SQUARE8;
}

void Square8Generator::pushElements(std::vector<esint> &elements, const std::vector<esint> &indices) const
{
    elements.push_back(indices[0]);
    elements.push_back(indices[2]);
    elements.push_back(indices[8]);
    elements.push_back(indices[6]);

    elements.push_back(indices[1]);
    elements.push_back(indices[5]);
    elements.push_back(indices[7]);
    elements.push_back(indices[3]);
}



