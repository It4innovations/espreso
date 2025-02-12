
#include "nglib.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.h"
#include "wrappers/nglib/w.nglib.h"

using namespace espreso;

NGLibGenerator::NGLibGenerator(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void NGLibGenerator::load()
{
    eslog::startln("NGLIB GENERATOR: STARTED", "NGLIB GENERATOR");
    profiler::syncstart("nglib");

    NGLib::generate(*this);

    body.resize(etype.size());
    material.resize(etype.size());
    profiler::syncend("nglib");
    eslog::endln("NGLIB GENERATOR: NGLIB GENERATED");
}



