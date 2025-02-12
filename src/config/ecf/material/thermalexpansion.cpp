
#include "thermalexpansion.h"
#include "config/configuration.hpp"

espreso::ThermalExpansionConfiguration::ThermalExpansionConfiguration()
: model(MODEL::ISOTROPIC),
  thermal_expansion(3)
{
    REGISTER(model, ECFMetaData()
            .setdescription({ "Tensors model." })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("ISOTROPIC").setdescription("Isotropic model."))
            .addoption(ECFOption().setname("ORTHOTROPIC").setdescription("Orthotropic model.")));

    ecfdescription->addSeparator();

    ecfdescription->registerParameter("TEX", thermal_expansion.get(0, 0), ECFMetaData()
            .setdescription({ "Thermal expansion X." })
            .setdatatype({ ECFDataType::EXPRESSION })
            .settensor(thermal_expansion));

    ecfdescription->registerParameter("TEY", thermal_expansion.get(1, 1), ECFMetaData()
            .setdescription({ "Thermal expansion Y." })
            .setdatatype({ ECFDataType::EXPRESSION })
            .settensor(thermal_expansion)
            .allowonly([&] () { return model != MODEL::ORTHOTROPIC; }));
    ecfdescription->registerParameter("TEZ", thermal_expansion.get(2, 2), ECFMetaData()
            .setdescription({ "Thermal expansion Z." })
            .setdatatype({ ECFDataType::EXPRESSION })
            .settensor(thermal_expansion)
            .allowonly([&] () { return model != MODEL::ORTHOTROPIC; }));
}

