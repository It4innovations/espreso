
#include "cm.h"

#include "esel.h"
#include "nsel.h"

#include "input/parsers/ansyscdb/ansyscdb.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include <cstring>
#include <functional>
#include <algorithm>

using namespace espreso;

CM::CM()
: entity(Entity::ELEMENTS)
{
    memset(name, '\0', MAX_NAME_SIZE);
}

CM& CM::parse(const char* begin)
{
    std::string commandLine = Parser::removecomments(Parser::getLine(begin), "!");
    std::vector<std::string> command = Parser::split(Parser::strip(commandLine), ",", false);

    switch (command.size()) {
    case 3:
        if (StringCompare::caseInsensitiveEq("VOLU", command[2])) {
            entity = Entity::VOLUME;
            eslog::error("ESPRESO Workbench parser error: implement cm, Entity='%s'\n", command[2].c_str());
        }
        if (StringCompare::caseInsensitiveEq("AREA", command[2])) {
            entity = Entity::AREA;
            eslog::error("ESPRESO Workbench parser error: implement cm, Entity='%s'\n", command[2].c_str());
        }
        if (StringCompare::caseInsensitiveEq("LINE", command[2])) {
            entity = Entity::LINE;
            eslog::error("ESPRESO Workbench parser error: implement cm, Entity='%s'\n", command[2].c_str());
        }
        if (StringCompare::caseInsensitiveEq("KP", command[2])) {
            entity = Entity::KP;
            eslog::error("ESPRESO Workbench parser error: implement cm, Entity='%s'\n", command[2].c_str());
        }
        if (StringCompare::caseInsensitiveEq("ELEM", command[2])) {
            entity = Entity::ELEMENTS;
        }
        if (StringCompare::caseInsensitiveEq("NODE", command[2])) {
            entity = Entity::NODES;
            eslog::error("ESPRESO Workbench parser error: implement cm, Entity='%s'\n", command[2].c_str());
        }
    case 2:
        memcpy(name, command[1].data(), command[1].size() < MAX_NAME_SIZE ? command[1].size() : MAX_NAME_SIZE);
        break;
    default:
        eslog::error("ESPRESO Workbench parser error: unknown format of '%s'\n", commandLine.c_str());
    }

    WorkbenchParser::fillIndices(begin, begin, begin);
    return *this;
}

bool CM::addRegion(
        const AnsysCDBData &mesh,
        const std::vector<ESel> &esel, std::map<std::string, std::vector<esint> > &eregions,
        const std::vector<NSel> &nsel, std::map<std::string, std::vector<esint> > &nregions)
{
    switch (entity) {
    case Entity::ELEMENTS:
        return addElementRegion(mesh, esel, eregions);
    case Entity::NODES:
        return addNodeRegion(nsel, nregions);
    default:
        return false;
    }
}

bool CM::addElementRegion(const AnsysCDBData &mesh, const std::vector<ESel> &esel, std::map<std::string, std::vector<esint> > &eregions)
{
    // clear relevant set means all elements
    std::vector<ESel> relevant;
    for (size_t i = 0; i < esel.size() && esel[i].first < first; i++) {
        switch (esel[i].type) {
        case ESel::Type::S:
        case ESel::Type::NONE:
            relevant.clear();
            relevant.push_back(esel[i]);
            break;
        case ESel::Type::A:
            if (relevant.size() == 1 && relevant[0].type == ESel::Type::NONE) {
                relevant.clear();
            }
            relevant.push_back(esel[i]);
            break;
        case ESel::Type::ALL:
            relevant.clear();
            break;
        default:
            if (fRank == info::mpi::rank) {
                eslog::error("ESPRESO Workbench parser error: not implemented type of esel command '%s'\n", esel[i].command().c_str());
            }
        }
    }

    if (relevant.size() == 0) {
        return true;
    }

    size_t threads = info::env::OMP_NUM_THREADS;
    std::vector<esint> edistribution = tarray<esint>::distribute(threads, mesh.esize.size());
    std::vector<std::vector<esint> > eid(threads);

    auto checkElements = [&] (std::function<void(size_t, esint)> chck) {
        #pragma omp parallel for
        for (size_t t = 0; t < threads; t++) {
            for (esint e = edistribution[t]; e < edistribution[t + 1]; ++e) {
                chck(t, e);
            }
        }
    };

    std::vector<esint> &data = eregions[name];
    for (size_t i = 0; i < relevant.size(); i++) {
        switch (relevant[i].item) {

        case ESel::Item::TYPE:
            if (relevant[i].VMIN == relevant[i].VMAX) {
                checkElements([&] (size_t t, esint e) {
                    if (mesh.et[e] == relevant[i].VMIN) {
                        eid[t].push_back(mesh.eIDs[e]);
                    }
                });
                break;
            }
            if (relevant[i].VINC == 1) {
                checkElements([&] (size_t t, esint e) {
                    if (relevant[i].VMIN <= mesh.et[e] && mesh.et[e] <= relevant[i].VMAX) {
                        eid[t].push_back(mesh.eIDs[e]);
                    }
                });
                break;
            }
            checkElements([&] (size_t t, esint e) {
                if (relevant[i].VMIN <= mesh.et[e] && mesh.et[e] <= relevant[i].VMAX && (mesh.et[e] - relevant[i].VMIN) % relevant[i].VINC == 0) {
                    eid[t].push_back(mesh.eIDs[e]);
                }
            });
            break;

        case ESel::Item::MAT:
            if (relevant[i].VMIN == relevant[i].VMAX) {
                checkElements([&] (size_t t, esint e) {
                    if (mesh.material[e] == relevant[i].VMIN) {
                        eid[t].push_back(mesh.eIDs[e]);
                    }
                });
                break;
            }
            if (relevant[i].VINC == 1) {
                checkElements([&] (size_t t, esint e) {
                    if (relevant[i].VMIN <= mesh.material[e] && mesh.material[e] <= relevant[i].VMAX) {
                        eid[t].push_back(mesh.eIDs[e]);
                    }
                });
                break;
            }
            checkElements([&] (size_t t, esint e) {
                if (relevant[i].VMIN <= mesh.material[e] && mesh.material[e] <= relevant[i].VMAX && (mesh.material[e] - relevant[i].VMIN) % relevant[i].VINC == 0) {
                    eid[t].push_back(mesh.eIDs[e]);
                }
            });
            break;

        default:
            if (fRank == info::mpi::rank) {
                eslog::error("ESPRESO Workbench parser error: not implemented type of esel command '%s'\n", esel[i].command().c_str());
            }
        }
    }

    for (size_t t = 0; t < threads; t++) {
        data.insert(data.end(), eid[t].begin(), eid[t].end());
    }
    utils::sortAndRemoveDuplicates(data);

    return true;
}

bool CM::addNodeRegion(const std::vector<NSel> &nsel, std::map<std::string, std::vector<esint> > &nregions)
{
    return false;
}
