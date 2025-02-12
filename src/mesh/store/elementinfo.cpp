
#include "elementinfo.h"

#include <algorithm>

using namespace espreso;

ElementsDistributionInfo::ElementsDistributionInfo()
: threads({ 0, 0 }), code(static_cast<int>(Element::CODE::SIZE))
{

}

void ElementsDistributionInfo::clear()
{
    process = DistributedDataInfo();
    std::fill(code.begin(), code.end(), DistributedDataInfo());
}
