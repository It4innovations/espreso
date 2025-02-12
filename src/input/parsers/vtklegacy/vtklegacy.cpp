
#include "vtklegacy.h"
#include "parser/geometry.h"

#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/parser.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "input/meshbuilder.h"

using namespace espreso;

VTKLegacyLoader::VTKLegacyLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void VTKLegacyLoader::load()
{
    eslog::startln("VTK PARSER: STARTED", "VTK PARSER");
    profiler::syncstart("vtk_legacy");

    std::vector<std::string> filepaths;
    std::vector<std::string> names;
    if (StringCompare::contains(_configuration.path, { "*" })) {
        size_t star = _configuration.path.find_last_of("*");
        std::string prefix = _configuration.path.substr(0, star);
        std::string suffix = _configuration.path.substr(star + 1);

        std::vector<std::string> allfiles;
        std::vector<char> serialized;
        std::string dir = utils::getFileDirectory(_configuration.path);
        if (info::mpi::rank == 0) {
            utils::listDirectory(dir, allfiles);
            serialized.push_back(0);
            for (auto file = allfiles.begin(); file != allfiles.end(); ++file) {
                ++serialized[0];
                serialized.push_back(file->size());
                serialized.insert(serialized.end(), file->data(), file->data() + file->size());
            }
        }
        Communication::broadcastUnknownSize(serialized);
        for (int i = 0, j = 1; i < serialized[0]; ++i, j += serialized[j] + 1) {
            std::string file = (dir[0] == '.' ?  "" : dir + "/") + std::string(serialized.data() + j + 1, serialized.data() + j + 1 + serialized[j]);
            if (StringCompare::caseSensitivePreffix(prefix, file) && StringCompare::caseSensitiveSuffix(file, suffix)) {
                filepaths.push_back(file);
                names.push_back(file.substr(prefix.size(), file.size() - prefix.size() - suffix.size()));
            }
        }
    } else {
        size_t nameend = _configuration.path.find_last_of(".");
        size_t namebegin = _configuration.path.find_last_of(".", nameend - 1) + 1;
        std::string name = _configuration.path.substr(namebegin, nameend - namebegin);
        filepaths.push_back(_configuration.path);
        names.push_back(name);
    }
    profiler::synccheckpoint("list_files");
    eslog::checkpointln("VTK PARSER: FILES LISTED");

    if (filepaths.size() == 0) {
        eslog::globalerror("VTK PARSER: incorrect path to VTK files: '%s'\n", _configuration.path.c_str());
    }

    InputFilePack pack;
    pack.commitFiles(filepaths);
    pack.prepare();
    profiler::synccheckpoint("prepare_reader");
    eslog::checkpointln("VTK PARSER: GEOMETRY READER PREPARED");

    pack.read();
    profiler::synccheckpoint("read");
    eslog::checkpointln("VTK PARSER: GEOMETRY READ");

    VTKLegacyGeometry geometry(pack);
    geometry.scan();
    profiler::synccheckpoint("scan");
    eslog::checkpointln("VTK PARSER: GEOMETRY SCANNED");

    geometry.parse(*this, names);
    removeDuplicates = true;
    body.resize(etype.size());
    material.resize(etype.size());
    profiler::synccheckpoint("parse");
    pack.clear();
    profiler::syncend("vtk_legacy");
    eslog::endln("VTK PARSER: GEOMETRY PARSED");
}
