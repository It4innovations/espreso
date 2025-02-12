
#include "geometry.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "mesh/element.h"
#include "input/meshbuilder.h"
#include "input/parsers/asciiparser.h"
#include "input/parsers/distributedscanner.h"

#include <numeric>

using namespace espreso;

EnsightGeometry::EnsightGeometry(InputFilePack &geofile)
: _geofile(geofile), _header{ Format::UNKNOWN, IDs::UNKNOWN, IDs::UNKNOWN }
{

}

void EnsightGeometry::scan()
{
    header();

    if (_header.format == Format::ASCII) {
        scanASCII();
    } else {
        scanBinary();
    }
}

void EnsightGeometry::parse(MeshBuilder &mesh)
{
    if (_header.format == Format::ASCII) {
        parseASCII(mesh);
    } else {
        parseBinary(mesh);
    }
}

void EnsightGeometry::header()
{
    // root always has header
    if (info::mpi::rank == 0) {
        const char *ndesc = _geofile.begin, *edesc = _geofile.begin;
        if (memcmp(_geofile.begin, "C Binary", 8) == 0) {
            _header.format = Format::BINARY;
            ndesc += 240, edesc += 320;
        } else {
            _header.format = Format::ASCII;
            while (*ndesc++ != '\n');
            while (*ndesc++ != '\n');
            edesc = ndesc;
            while (*edesc++ != '\n');
        }
        if (DistributedScanner::check(ndesc, "node id ")) {
            if (DistributedScanner::check(ndesc + 8, "off")) { _header.nodeIDs = IDs::OFF; }
            if (DistributedScanner::check(ndesc + 8, "given")) { _header.nodeIDs = IDs::GIVEN; }
            if (DistributedScanner::check(ndesc + 8, "assign")) { _header.nodeIDs = IDs::ASSIGN; }
            if (DistributedScanner::check(ndesc + 8, "ignore")) { _header.nodeIDs = IDs::INGNORE; }
        }
        if (DistributedScanner::check(edesc, "element id ")) {
            if (DistributedScanner::check(edesc + 11, "off")) { _header.elementIDs = IDs::OFF; }
            if (DistributedScanner::check(edesc + 11, "given")) { _header.elementIDs = IDs::GIVEN; }
            if (DistributedScanner::check(edesc + 11, "assign")) { _header.elementIDs = IDs::ASSIGN; }
            if (DistributedScanner::check(edesc + 11, "ignore")) { _header.elementIDs = IDs::INGNORE; }
        }
    }

    Communication::broadcast(&_header, sizeof(Header), MPI_BYTE, 0);

    if (_header.format == Format::UNKNOWN || _header.nodeIDs == IDs::UNKNOWN || _header.elementIDs == IDs::UNKNOWN) {
        eslog::globalerror("MESIO internal error: cannot recognize EnSight Gold geometry format.\n");
    }
}

void EnsightGeometry::scanBinary()
{
    DistributedScanner parser;

    auto getsize = [] (const void *data) {
        int size;
        memcpy(&size, data, sizeof(int));
        return size;
    };

    auto addcoordinates = [&] (const char *c) {
        _coordinates.push_back(Coordinates(_geofile.distribution[info::mpi::rank] + (c - _geofile.begin) + 80 + sizeof(int), getsize(c + 80)));
    };
    auto skipcoordinates = [&] (const char *c) -> size_t {
        int nn; memcpy(&nn, c + 80, sizeof(int));
        return 3 * nn * sizeof(float);
    };

    auto addelement = [&] (const char *c, Elements::Type type) {
        size_t desc = 80 + sizeof(int), offset = _geofile.distribution[info::mpi::rank] + (c - _geofile.begin);
        int nn; memcpy(&nn, c + 80, sizeof(int));
        _elements.push_back(Elements(type, offset + desc, nn));
    };

    auto skipelements = [&] (const char *c, int enodes) -> size_t {
        int nn; memcpy(&nn, c + 80, sizeof(int));
        return enodes * nn * sizeof(int);
    };

    parser.add("part"       , [&] (const char *c) { _parts.insert(_parts.end(), c + 80 + sizeof(int), c + 80 + sizeof(int) + 80); });
    parser.add("coordinates", [&] (const char *c) { addcoordinates(c); }                       , [&] (const char *c) { return skipcoordinates(c); });

    parser.add("point"      , [&] (const char *c) { addelement(c, Elements::Type::POINT); }    , [&] (const char *c) { return skipelements(c,  1); });
    parser.add("bar2"       , [&] (const char *c) { addelement(c, Elements::Type::BAR2); }     , [&] (const char *c) { return skipelements(c,  2); });
    parser.add("bar3"       , [&] (const char *c) { addelement(c, Elements::Type::BAR3); }     , [&] (const char *c) { return skipelements(c,  3); });
    parser.add("tria3"      , [&] (const char *c) { addelement(c, Elements::Type::TRIA3); }    , [&] (const char *c) { return skipelements(c,  3); });
    parser.add("tria6"      , [&] (const char *c) { addelement(c, Elements::Type::TRIA6); }    , [&] (const char *c) { return skipelements(c,  6); });
    parser.add("quad4"      , [&] (const char *c) { addelement(c, Elements::Type::QUAD4); }    , [&] (const char *c) { return skipelements(c,  4); });
    parser.add("quad8"      , [&] (const char *c) { addelement(c, Elements::Type::QUAD8); }    , [&] (const char *c) { return skipelements(c,  8); });
    parser.add("tetra4"     , [&] (const char *c) { addelement(c, Elements::Type::TETRA4); }   , [&] (const char *c) { return skipelements(c,  4); });
    parser.add("tetra10"    , [&] (const char *c) { addelement(c, Elements::Type::TETRA10); }  , [&] (const char *c) { return skipelements(c, 10); });
    parser.add("pyramid5"   , [&] (const char *c) { addelement(c, Elements::Type::PYRAMID5); } , [&] (const char *c) { return skipelements(c,  5); });
    parser.add("pyramid13"  , [&] (const char *c) { addelement(c, Elements::Type::PYRAMID13); }, [&] (const char *c) { return skipelements(c, 13); });
    parser.add("penta6"     , [&] (const char *c) { addelement(c, Elements::Type::PENTA6); }   , [&] (const char *c) { return skipelements(c,  6); });
    parser.add("penta15"    , [&] (const char *c) { addelement(c, Elements::Type::PENTA15); }  , [&] (const char *c) { return skipelements(c, 15); });
    parser.add("hexa8"      , [&] (const char *c) { addelement(c, Elements::Type::HEXA8); }    , [&] (const char *c) { return skipelements(c,  8); });
    parser.add("hexa20"     , [&] (const char *c) { addelement(c, Elements::Type::HEXA20); }   , [&] (const char *c) { return skipelements(c, 20); });

    parser.add("nsided"     , [&] (const char *c) { eslog::error("Ensight parser error: ESPRESO does not support nsided elements.\n"); });
    parser.add("nfaced"     , [&] (const char *c) { eslog::error("Ensight parser error: ESPRESO does not support nfaced elements.\n"); });

    parser.scan(_geofile);
    parser.synchronize(_parts, _coordinates, _elements);
}

void EnsightGeometry::parseBinary(MeshBuilder &mesh)
{
    // skip ids: currently we always generate them
    if (_header.nodeIDs == IDs::GIVEN || _header.nodeIDs == IDs::INGNORE) {
        for (size_t i = 0; i < _coordinates.size(); ++i) {
            _coordinates[i].offset += _coordinates[i].nn * sizeof(int);
        }
    }
    if (_header.elementIDs == IDs::GIVEN || _header.elementIDs == IDs::INGNORE) {
        for (size_t i = 0; i < _elements.size(); ++i) {
            _elements[i].offset += _elements[i].ne * sizeof(int);
        }
    }

    auto round = [] (size_t &current, size_t &start, size_t chunk) {
        if ((current - start) % chunk != 0) {
            current += chunk - (current - start) % chunk;
        }
    };

    struct eoffsetinfo {
        std::string name;
        int offset;
        int size;
        int totalsize;
    };

    std::vector<eoffsetinfo> einfo;
    size_t parts = _parts.size() / 80;
    esint cidoffset = 0;

    for (size_t p = 0, i = 0; p < parts; ++p) {
        auto begin = _parts.begin() + p * 80;
        auto end = _parts.begin() + p * 80;
        while (*end != 0) { ++end; }
        std::string name = std::string(begin, end);

        int firstrank = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[p].offset) - _geofile.distribution.begin() - 1;
        int lastrank = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[p].offset + 3 * sizeof(float) * _coordinates[p].nn) - _geofile.distribution.begin() - 1;

        size_t cbegin = std::max(_coordinates[p].offset, _geofile.distribution[info::mpi::rank]);
        size_t cend = std::min(_coordinates[p].offset + 3 * sizeof(float) * _coordinates[p].nn, _geofile.distribution[info::mpi::rank + 1]);
        if (cbegin < cend) {
            // read data: x x x x x...; y y y y y y...; z z z z z z...
            round(cbegin, _coordinates[p].offset, sizeof(float));
            round(cend, _coordinates[p].offset, sizeof(float));
            std::vector<float> coordinates((cend - cbegin) / sizeof(float));
            memcpy(coordinates.data(), _geofile.begin + cbegin - _geofile.distribution[info::mpi::rank], cend - cbegin);

            std::vector<size_t> cdistribution = tarray<size_t>::distribute(lastrank - firstrank + 1, _coordinates[p].nn);
            if (firstrank != lastrank) { // group data: x y z; x y z; x y z; ...
                cbegin = (cbegin - _coordinates[p].offset) / sizeof(float); // from char to float offset
                cend   = (cend   - _coordinates[p].offset) / sizeof(float); // from char to float offset
                std::vector<int> sBuffer, rBuffer;
                for (int r = firstrank; r <= lastrank; ++r) {
                    size_t prevsize = sBuffer.size();
                    sBuffer.push_back(0); // total size
                    sBuffer.push_back(r); // target
                    sBuffer.push_back(info::mpi::rank); // source
                    sBuffer.push_back(0); // x size
                    sBuffer.push_back(0); // y size
                    sBuffer.push_back(0); // z size

                    for (int d = 0; d < 3; ++d) {
                        size_t rbegin = std::max(d * _coordinates[p].nn + cdistribution[r - firstrank], cbegin);
                        size_t rend   = std::min(d * _coordinates[p].nn + cdistribution[r - firstrank + 1], cend);
                        if (rend < rbegin) {
                            rend = rbegin;
                        }
                        sBuffer[prevsize + 3 + d] = rend - rbegin;
                        if (rbegin != rend) {
                            sBuffer.insert(sBuffer.end(), reinterpret_cast<const int*>(coordinates.data() + rbegin - cbegin), reinterpret_cast<const int*>(coordinates.data() + rend - cbegin));
                        }
                    }
                    sBuffer[prevsize] = sBuffer.size() - prevsize;
                }

                if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer, firstrank, lastrank + 1)) {
                    eslog::error("Ensight Gold parser: cannot exchange coordinate data.\n");
                }

                coordinates.resize(3 * (cdistribution[info::mpi::rank - firstrank + 1] - cdistribution[info::mpi::rank - firstrank]));
                std::vector<int> alloffsets(3 * (lastrank - firstrank + 1)); // we have to keep data order
                size_t offset = 0;
                for (int r = firstrank; r <= lastrank; ++r) {
                    ++offset; // totalsize
                    ++offset; // target (always me)
                    int source = rBuffer[offset++] - firstrank;
                    for (int d = 0; d < 3; ++d) {
                        alloffsets[3 * source + d] = rBuffer[offset++];
                    }
                    for (int d = 0; d < 3; ++d) {
                        offset += alloffsets[3 * source + d];
                    }
                }
                for (int d = 0; d < 3; ++d) {
                    int sum = d * (cdistribution[info::mpi::rank - firstrank + 1] - cdistribution[info::mpi::rank - firstrank]);
                    for (int r = 0; r <= lastrank - firstrank; ++r) {
                        int tmp = alloffsets[3 * r + d];
                        alloffsets[3 * r + d] = sum;
                        sum += tmp;
                    }
                }
                offset = 0;
                for (int r = firstrank; r <= lastrank; ++r) {
                    ++offset; // totalsize
                    ++offset; // target (always me)
                    int source = rBuffer[offset++] - firstrank;
                    int xyzsize[3] = { rBuffer[offset++], rBuffer[offset++], rBuffer[offset++] };
                    for (int d = 0; d < 3; ++d) {
                        if (xyzsize[d]) {
                            memcpy(coordinates.data() + alloffsets[3 * source + d], rBuffer.data() + offset, sizeof(float) * xyzsize[d]);
                            offset += xyzsize[d] * (sizeof(float) / sizeof(int));
                        }
                    }
                }
            }

            int ncoordinates = coordinates.size() / 3;
            mesh.coordinates.reserve(mesh.coordinates.size() + ncoordinates);
            for (int c = 0; c < ncoordinates; ++c) {
                mesh.coordinates.push_back(Point(coordinates[ncoordinates * 0 + c], coordinates[ncoordinates * 1 + c], coordinates[ncoordinates * 2 + c]));
            }
            mesh.nIDs.resize(mesh.coordinates.size());
            std::iota(mesh.nIDs.end() - ncoordinates, mesh.nIDs.end(), cidoffset + cdistribution[info::mpi::rank - firstrank]);
        }

        int ntotalelements = 0;
        for (; i < _elements.size() && (p + 1 == parts || _elements[i].offset < _coordinates[p + 1].offset); ++i) {
            int esize; Element::CODE code;
            switch (_elements[i].type) {
            case Elements::Type::POINT    : esize =  1; code = Element::CODE::POINT1   ; break;
            case Elements::Type::BAR2     : esize =  2; code = Element::CODE::LINE2    ; break;
            case Elements::Type::BAR3     : esize =  3; code = Element::CODE::LINE3    ; break;
            case Elements::Type::TRIA3    : esize =  3; code = Element::CODE::TRIANGLE3; break;
            case Elements::Type::TRIA6    : esize =  6; code = Element::CODE::TRIANGLE6; break;
            case Elements::Type::QUAD4    : esize =  4; code = Element::CODE::SQUARE4  ; break;
            case Elements::Type::QUAD8    : esize =  8; code = Element::CODE::SQUARE8  ; break;
            case Elements::Type::TETRA4   : esize =  4; code = Element::CODE::TETRA4   ; break;
            case Elements::Type::TETRA10  : esize = 10; code = Element::CODE::TETRA10  ; break;
            case Elements::Type::PYRAMID5 : esize =  5; code = Element::CODE::PYRAMID5 ; break;
            case Elements::Type::PYRAMID13: esize = 13; code = Element::CODE::PYRAMID13; break;
            case Elements::Type::PENTA6   : esize =  6; code = Element::CODE::PRISMA6  ; break;
            case Elements::Type::PENTA15  : esize = 15; code = Element::CODE::PRISMA15 ; break;
            case Elements::Type::HEXA8    : esize =  8; code = Element::CODE::HEXA8    ; break;
            case Elements::Type::HEXA20   : esize = 20; code = Element::CODE::HEXA20   ; break;
            case Elements::Type::NSIDED   : esize =  0; code = Element::CODE::SIZE     ; break; // not supported
            case Elements::Type::NFACED   : esize =  0; code = Element::CODE::SIZE     ; break; // not supported
            default: esize = 0; code = Element::CODE::SIZE;
            }

            size_t ebegin = std::max(_elements[i].offset, _geofile.distribution[info::mpi::rank]);
            size_t eend = std::min(_elements[i].offset + esize * sizeof(int) * _elements[i].ne, _geofile.distribution[info::mpi::rank + 1]);
            if (ebegin < eend) {
                round(ebegin, _elements[i].offset, esize * sizeof(int));
                round(eend, _elements[i].offset, esize * sizeof(int));
                std::vector<int> elements((eend - ebegin) / sizeof(int));
                memcpy(elements.data(), _geofile.begin + ebegin - _geofile.distribution[info::mpi::rank], eend - ebegin);

                int nelements = (eend - ebegin) / (sizeof(int) * esize);
                if (code != Element::CODE::POINT1) {
                    ntotalelements += nelements;
                    mesh.etype.resize(mesh.etype.size() + nelements, (int)code);
                    mesh.esize.resize(mesh.esize.size() + nelements, esize);
                    mesh.enodes.reserve(mesh.enodes.size() + elements.size());
                    for (size_t n = 0; n < elements.size(); ++n) {
                        mesh.enodes.push_back(elements[n] + cidoffset - 1);
                    }
                } else {
                    for (size_t n = 0; n < elements.size(); ++n) {
                        elements[n] += cidoffset - 1;
                    }
                    mesh.nregions[name].assign(elements.begin(), elements.end());
                }
            } else {
                if (code == Element::CODE::POINT1) {
                    mesh.nregions[name] = {};
                }
            }
        }
        einfo.push_back(eoffsetinfo{.name=name, .offset=ntotalelements, .size=ntotalelements, .totalsize=0});
        cidoffset += _coordinates[p].nn;
    }

    mesh.eIDs.resize(mesh.etype.size());

    std::vector<esint> sum, offset;
    for (size_t i = 0; i < einfo.size(); ++i) {
        offset.push_back(einfo[i].offset);
    }
    sum.resize(offset.size());
    Communication::exscan(sum, offset);
    for (size_t i = 0; i < einfo.size(); ++i) {
        einfo[i].offset = offset[i];
        einfo[i].totalsize = sum[i];
    }
    for (size_t i = 0, prevsize = 0, eidoffset = 0; i < einfo.size(); ++i) {
        if (einfo[i].totalsize) { // non-node region
            std::iota(mesh.eIDs.begin() + prevsize, mesh.eIDs.begin() + prevsize + einfo[i].size, einfo[i].offset + eidoffset);
            mesh.eregions[einfo[i].name] = std::vector<esint>(mesh.eIDs.begin() + prevsize, mesh.eIDs.begin() + prevsize + einfo[i].size);
            prevsize += einfo[i].size;
            eidoffset += einfo[i].totalsize;
        }
    }
}

void EnsightGeometry::scanASCII()
{
    DistributedScanner parser;

    auto addelement = [&] (const char *c, Elements::Type type) {
        while (*c++ != '\n');
        int nn = atof(c);
        while (*c++ != '\n');
        int offset = _geofile.distribution[info::mpi::rank] + (c - _geofile.begin);
        _elements.push_back(Elements(type, offset, nn));
    };

    auto skipelements = [&] (const char *c, int enodes) -> int {
        while (*c++ != '\n');
        int nn = atof(c);
        return (enodes * 10 + 1) * nn;
    };

    auto addcoordinates = [&] (const char *c) {
        while (*c++ != '\n');
        int count = atof(c);
        while (*c++ != '\n');
        _coordinates.push_back(Coordinates(_geofile.distribution[info::mpi::rank] + (c - _geofile.begin), count));
    };

    auto skipcoordinates = [&] (const char *c) -> size_t {
        while (*c++ != '\n');
        int count = atof(c);
        while (*c++ != '\n');
        return 13 * count;
    };

    auto addparts = [&] (const char *c, std::vector<char> &vector) {
        while (*c++ != '\n');
        while (*c++ != '\n');
        size_t size = vector.size();
        vector.resize(size + 80);
        while (*c != '\n'){
            vector[size++] = *c++;
        }
    };

    parser.add("part"       , [&] (const char *c) { addparts(c, _parts); });
    parser.add("coordinates", [&] (const char *c) { addcoordinates(c); }                       , [&] (const char *c) { return skipcoordinates(c); });

    parser.add("point"      , [&] (const char *c) { addelement(c, Elements::Type::POINT); }    , [&] (const char *c) { return skipelements(c,  1); });
    parser.add("bar2"       , [&] (const char *c) { addelement(c, Elements::Type::BAR2); }     , [&] (const char *c) { return skipelements(c,  2); });
    parser.add("bar3"       , [&] (const char *c) { addelement(c, Elements::Type::BAR3); }     , [&] (const char *c) { return skipelements(c,  3); });
    parser.add("tria3"      , [&] (const char *c) { addelement(c, Elements::Type::TRIA3); }    , [&] (const char *c) { return skipelements(c,  3); });
    parser.add("tria6"      , [&] (const char *c) { addelement(c, Elements::Type::TRIA6); }    , [&] (const char *c) { return skipelements(c,  6); });
    parser.add("quad4"      , [&] (const char *c) { addelement(c, Elements::Type::QUAD4); }    , [&] (const char *c) { return skipelements(c,  4); });
    parser.add("quad8"      , [&] (const char *c) { addelement(c, Elements::Type::QUAD8); }    , [&] (const char *c) { return skipelements(c,  8); });
    parser.add("tetra4"     , [&] (const char *c) { addelement(c, Elements::Type::TETRA4); }   , [&] (const char *c) { return skipelements(c,  4); });
    parser.add("tetra10"    , [&] (const char *c) { addelement(c, Elements::Type::TETRA10); }  , [&] (const char *c) { return skipelements(c,  10); });
    parser.add("pyramid5"   , [&] (const char *c) { addelement(c, Elements::Type::PYRAMID5); } , [&] (const char *c) { return skipelements(c,  5); });
    parser.add("pyramid13"  , [&] (const char *c) { addelement(c, Elements::Type::PYRAMID13); }, [&] (const char *c) { return skipelements(c, 13); });
    parser.add("penta6"     , [&] (const char *c) { addelement(c, Elements::Type::PENTA6); }   , [&] (const char *c) { return skipelements(c,  6); });
    parser.add("penta15"    , [&] (const char *c) { addelement(c, Elements::Type::PENTA15); }  , [&] (const char *c) { return skipelements(c, 15); });
    parser.add("hexa8"      , [&] (const char *c) { addelement(c, Elements::Type::HEXA8); }    , [&] (const char *c) { return skipelements(c,  8); });
    parser.add("hexa20"     , [&] (const char *c) { addelement(c, Elements::Type::HEXA20); }   , [&] (const char *c) { return skipelements(c, 20); });

    parser.add("nsided"     , [&] (const char *c) { eslog::error("Ensight parser error: ESPRESO does not support nsided elements.\n"); });
    parser.add("nfaced"     , [&] (const char *c) { eslog::error("Ensight parser error: ESPRESO does not support nfaced elements.\n"); });

    parser.scan(_geofile);
    parser.synchronize(_parts, _coordinates, _elements);
}

void EnsightGeometry::parseASCII(MeshBuilder &mesh)
{
    // skip ids: currently we always generate them
    if (_header.nodeIDs == IDs::GIVEN || _header.nodeIDs == IDs::INGNORE) {
        for (size_t i = 0; i < _coordinates.size(); ++i) {
            _coordinates[i].offset += _coordinates[i].nn * 11;
        }
    }
    if (_header.elementIDs == IDs::GIVEN || _header.elementIDs == IDs::INGNORE) {
        for (size_t i = 0; i < _elements.size(); ++i) {
            _elements[i].offset += _elements[i].ne * 11;
        }
    }

    auto align = [] (size_t current, size_t begin) {
        return (current - begin) % 13 == 0 ? current : current + 13 - ((current - begin) % 13);
    };

    for (size_t i = 0, offset = 0; i < _coordinates.size(); offset += _coordinates[i++].nn) {
        std::vector<float> coordinates;
        size_t A = std::max(_coordinates[i].offset, _geofile.distribution[info::mpi::rank]);
        size_t B = std::min(_coordinates[i].offset + 3 * 13 * _coordinates[i].nn, _geofile.distribution[info::mpi::rank + 1]);
        if (A < B) {
            A = align(A, _coordinates[i].offset) - _geofile.distribution[info::mpi::rank];
            B = align(B, _coordinates[i].offset) - _geofile.distribution[info::mpi::rank];
            for (const char *c = _geofile.begin + A; c < _geofile.begin + B; c += 13) {
                coordinates.push_back(atof(c));
            }

            int start = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[i].offset) - _geofile.distribution.begin() - 1;
            int end = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[i].offset + 3 * 13 * _coordinates[i].nn) - _geofile.distribution.begin() - 1;

            if (start == info::mpi::rank) {
                std::vector<float, initless_allocator<float> > rbuffer(3 * _coordinates[i].nn);
                std::copy(coordinates.begin(), coordinates.end(), rbuffer.begin());
                std::vector<MPI_Request> req(end - start);
                for (int r = start + 1; r <= end; ++r) {
                    size_t rbegin = align(std::max(_coordinates[i].offset, _geofile.distribution[r]), _coordinates[i].offset) - _coordinates[i].offset;
                    size_t rend = align(std::min(_coordinates[i].offset + 3 * 13 * _coordinates[i].nn, _geofile.distribution[r + 1]), _coordinates[i].offset) - _coordinates[i].offset;
                    MPI_Irecv(rbuffer.data() + rbegin / 13, (rend - rbegin) / 13, MPI_FLOAT, r, 0, info::mpi::comm, req.data() + (r - start - 1));
                }
                mesh.nIDs.reserve(mesh.coordinates.size() + _coordinates[i].nn);
                mesh.coordinates.reserve(mesh.coordinates.size() + _coordinates[i].nn);
                MPI_Waitall(end - start, req.data(), MPI_STATUSES_IGNORE);
                for (int n = 0; n < _coordinates[i].nn; ++n) {
                    mesh.nIDs.push_back(offset + n);
                    mesh.coordinates.push_back(Point(rbuffer[n + _coordinates[i].nn * 0], rbuffer[n + _coordinates[i].nn * 1], rbuffer[n + _coordinates[i].nn * 2]));
                }
            } else {
                MPI_Send(coordinates.data(), coordinates.size(), MPI_FLOAT, start, 0, info::mpi::comm);
            }
        }
    }

    auto alignElement = [] (size_t current, size_t begin, int line) {
        return (current - begin) % line == 0 ? current : current + line - ((current - begin) % line);
    };

    size_t parts = _parts.size() / 80;
    size_t coordinate_id_offset = 0;
    size_t element_id_offset = 0;
    for (size_t p = 0, e = 0; p < parts; ++p) {
        auto begin = _parts.begin() + p * 80;
        auto end = _parts.begin() + p * 80;
        while (*end != 0) { ++end; }
        std::string name = std::string(begin, end);
        int esize;
        Element::CODE code = Element::CODE::SIZE;
        std::vector<esint> eids;
        size_t prevsize = mesh.eIDs.size();
        for (; e < _elements.size() && (p + 1 == parts || _elements[e].offset < _coordinates[p + 1].offset); ++e) {
            switch (_elements[e].type){
            case Elements::Type::POINT     : esize =  1; code = Element::CODE::POINT1;    break;
            case Elements::Type::BAR2      : esize =  2; code = Element::CODE::LINE2;     break;
            case Elements::Type::TRIA3     : esize =  3; code = Element::CODE::TRIANGLE3; break;
            case Elements::Type::QUAD4     : esize =  4; code = Element::CODE::SQUARE4;   break;
            case Elements::Type::TETRA4    : esize =  4; code = Element::CODE::TETRA4;    break;
            case Elements::Type::PYRAMID5  : esize =  5; code = Element::CODE::PYRAMID5;  break;
            case Elements::Type::PENTA6    : esize =  6; code = Element::CODE::PRISMA6;   break;
            case Elements::Type::HEXA8     : esize =  8; code = Element::CODE::HEXA8;     break;
            case Elements::Type::BAR3      : esize =  3; code = Element::CODE::LINE3;     break;
            case Elements::Type::TRIA6     : esize =  6; code = Element::CODE::TRIANGLE6; break;
            case Elements::Type::QUAD8     : esize =  8; code = Element::CODE::SQUARE8;   break;
            case Elements::Type::TETRA10   : esize = 10; code = Element::CODE::TETRA10;   break;
            case Elements::Type::PYRAMID13 : esize = 13; code = Element::CODE::PYRAMID13; break;
            case Elements::Type::PENTA15   : esize = 15; code = Element::CODE::PRISMA15;  break;
            case Elements::Type::HEXA20    : esize = 20; code = Element::CODE::HEXA20;    break;
            case Elements::Type::NSIDED    : esize =  0; code = Element::CODE::SIZE;      break; // not supported
            case Elements::Type::NFACED    : esize =  0; code = Element::CODE::SIZE;      break; // not supported
            default                        : esize =  0; code = Element::CODE::SIZE;
            }

            size_t E = std::max(_elements[e].offset, _geofile.distribution[info::mpi::rank]);
            size_t F = std::min(_elements[e].offset + _elements[e].ne * (esize * 10 + 1), _geofile.distribution[info::mpi::rank + 1]);
            if (E < F) {
                E = alignElement(E, _elements[e].offset, 10 * esize + 1);
                F = alignElement(F, _elements[e].offset, 10 * esize + 1);
                int eid_offset = (E - _elements[e].offset) / (10 * esize + 1);
                E -= _geofile.distribution[info::mpi::rank];
                F -= _geofile.distribution[info::mpi::rank];

                if (code == Element::CODE::POINT1){
                    for (const char *d = _geofile.begin + E; d < _geofile.begin + F; d += 11, ++eid_offset) {
                        eids.push_back(atof(d) - 1 + coordinate_id_offset);
                    }
                } else {
                    for (const char *d = _geofile.begin + E; d < _geofile.begin + F; d += 10 * esize + 1, ++eid_offset) {
                        mesh.eIDs.push_back(element_id_offset + eid_offset);
                        for (int n = 0; n < esize; ++n) {
                            mesh.enodes.push_back(atof(d + 10 * n) - 1 + coordinate_id_offset);
                        }
                    }
                    mesh.esize.resize(mesh.eIDs.size(), esize);
                    mesh.etype.resize(mesh.eIDs.size(), (int)code);
                }
            }
            element_id_offset += _elements[e].ne;
        }
        if (code == Element::CODE::POINT1){
            mesh.nregions[name] = eids;
        } else {
            mesh.eregions[name].assign(mesh.eIDs.begin() + prevsize, mesh.eIDs.end());
        }

        // go to another part
        coordinate_id_offset += _coordinates[p].nn;
    }
}
