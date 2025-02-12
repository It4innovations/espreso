
#include "xdmf.h"
#include "lightdata/lightdata.h"
#include "heavydata/griddata.h"
#include "basis/logging/profiler.h"
#include "esinfo/eslog.h"
#include "config/ecf/input/input.h"

using namespace espreso;


XDMFLoader::XDMFLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void XDMFLoader::load()
{
    eslog::startln("XDMF PARSER: STARTED", "XDMF PARSER");
    profiler::syncstart("xdmf");

    LightData lightdata(_configuration.path);
    profiler::synccheckpoint("lightdata");
    eslog::checkpointln("XDMF PARSER: LIGHTDATA PARSED");

    GridData grid(lightdata);
    grid.scan();
    profiler::synccheckpoint("scan");
    eslog::checkpointln("XDMF PARSER: LIGHTDATA SCANNED");

    grid.read();
    profiler::synccheckpoint("read");
    eslog::checkpointln("XDMF PARSER: HEAVYDATA READ");

    grid.parse(*this);
    removeDuplicates = true;
    body.resize(etype.size());
    material.resize(etype.size());
    profiler::synccheckpoint("parse");
    profiler::syncend("xdmf");
    eslog::endln("XDMF PARSER: HEAVYDATA PARSED");
}
