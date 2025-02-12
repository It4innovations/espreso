
#include "gridset.h"

#include "config/configuration.hpp"

espreso::GridSetGeneratorConfiguration::GridSetGeneratorConfiguration()
{
    REGISTER(grids, ECFMetaData()
            .setdescription({ "An index of grid in tower. Indices has to be continuous starting from 0.", "Description of grid." })
            .setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
            .setpattern({ "0" }));
}


