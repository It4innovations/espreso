
#include "dualgraph.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

void DualGraph::clear()
{
    vertices.clear();
    edges.clear();
}

void DualGraph::pushVertex(int index, int size)
{
    vertices[index].kernel.size = size;
}

void DualGraph::initVertices()
{
    int offset = 0;
    for (auto v = vertices.begin(); v != vertices.end(); ++v) {
        offset += v->second.kernel.size;
    }
    Communication::exscan(offset);
    for (auto v = vertices.begin(); v != vertices.end(); ++v) {
        v->second.rank = info::mpi::rank;
        v->second.kernel.goffset = offset;
        offset += v->second.kernel.size;
    }
}

void DualGraph::setFromDomains(const DecompositionFETI *decomposition, const std::vector<int> &lMap)
{
    for (size_t i = 0; i < lMap.size(); ) {
        int domains = lMap[i + 1];
        for (int d1 = 0; d1 < domains; ++d1) {
            for (int d2 = 0; d2 < domains; ++d2) {
                edges[lMap[i + 2 + d1]].push_back(lMap[i + 2 + d2]);
            }
        }
        i += lMap[i + 1] + 2;
    }
    for (auto e = edges.begin(); e != edges.end(); ++e) {
        utils::sortAndRemoveDuplicates(e->second);
    }

    std::vector<std::vector<int> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
    for (auto v = vertices.begin(); v != vertices.end(); ++v) {
        int last = -1;
        for (size_t ve = 0; ve < edges[v->first].size(); ++ve) {
            int n = decomposition->noffset(edges[v->first][ve]);
            if (!decomposition->ismy(edges[v->first][ve]) && last < n) {
                sBuffer[n].push_back(v->first);
                sBuffer[n].push_back(v->second.kernel.goffset);
                sBuffer[n].push_back(v->second.kernel.size);
                last = n;
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
        eslog::error("cannot exchange dual graph info\n");
    }

    for (size_t n = 0; n < rBuffer.size(); ++n) {
        for (size_t i = 0; i < rBuffer[n].size(); i += 3) {
            vertices[rBuffer[n][i]].rank = decomposition->neighbors[n];
            vertices[rBuffer[n][i]].kernel.goffset = rBuffer[n][i + 1];
            vertices[rBuffer[n][i]].kernel.size = rBuffer[n][i + 2];
        }
    }

    for (size_t i = 0, offset = 0; i < lMap.size(); ) {
        int lsize = lMap[i];
        int domains = lMap[i + 1];
        for (int d = 0; d < domains; ++d) {
            vertices[lMap[i + 2 + d]].lambdas.indices.push_back({ (int)offset, lsize });
            vertices[lMap[i + 2 + d]].lambdas.total += lsize;
            std::vector<int> &lneighs = vertices[lMap[i + 2 + d]].lambdas.indices.back().lower_neighs;
            if (decomposition->ismy(lMap[i + 2 + d])) {
                for (int lower = 0; lower < d; ++lower) {
                    if (lMap[i + 2 + lower] < decomposition->dbegin) {
                        int n = decomposition->noffset(lMap[i + 2 + lower]);
                        if (lneighs.empty() || lneighs.back() != n) {
                            lneighs.push_back(n);
                        }
                    }
                }
            }
        }
        i += lMap[i + 1] + 2;
        offset += lsize;
    }
}

void DualGraph::setFromClusters(const DecompositionFETI *decomposition, const std::vector<int> &lMap)
{

}

void DualGraph::spread(const DecompositionFETI *decomposition)
{
    std::map<int, std::vector<int> > next = edges;

    for (auto e = edges.begin(); e != edges.end(); ++e) {
        for (size_t ve = 0; ve < e->second.size(); ++ve) {
            for (size_t vt = 0; vt < edges[e->second[ve]].size(); ++vt) {
                next[e->first].push_back(edges[e->second[ve]][vt]);
            }
        }
    }

    for (auto e = next.begin(); e != next.end(); ++e) {
        utils::sortAndRemoveDuplicates(e->second);
    }

    std::vector<std::vector<int> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
    if (decomposition->neighbors.size()) {
        sBuffer[0].push_back(vertices.size());
        for (auto v = vertices.begin(); v != vertices.end(); ++v) {
            sBuffer[0].push_back(v->first);
            sBuffer[0].push_back(v->second.rank);
            sBuffer[0].push_back(v->second.kernel.goffset);
            sBuffer[0].push_back(v->second.kernel.size);
        }
        for (auto e = next.begin(); e != next.end(); ++e) {
            sBuffer[0].push_back(e->first);
            sBuffer[0].push_back(e->second.size());
            for (size_t ve = 0; ve < e->second.size(); ++ve) {
                sBuffer[0].push_back(e->second[ve]);
            }
        }
        for (size_t n = 1; n < sBuffer.size(); ++n) {
            sBuffer[n] = sBuffer[0];
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
        eslog::error("cannot spread domain dual graph.\n");
    }

    for (size_t n = 0; n < rBuffer.size(); ++n) {
        size_t i = 1;
        for (int v = 0; v < rBuffer[n][0]; ++v, i += 4) {
            vertices[rBuffer[n][i]].rank = rBuffer[n][i + 1];
            vertices[rBuffer[n][i]].kernel.goffset = rBuffer[n][i + 2];
            vertices[rBuffer[n][i]].kernel.size    = rBuffer[n][i + 3];
        }
        while (i < rBuffer[n].size()) {
            for (int e = 0; e < rBuffer[n][i + 1]; ++e) {
                next[rBuffer[n][i]].push_back(rBuffer[n][i + 2 + e]);
            }
            i += rBuffer[n][i + 1] + 2;
        }
    }

    for (auto e = next.begin(); e != next.end(); ++e) {
        utils::sortAndRemoveDuplicates(e->second);
    }
    edges.swap(next);
}

void DualGraph::print()
{
    Communication::serialize([&] () {
        printf(" --- rank  %2d ---\n", info::mpi::rank);
        for (auto v = vertices.begin(); v != vertices.end(); ++v) {
            printf("%2d :: holder:%2d, goffset:%2d\n", v->first, v->second.rank, v->second.kernel.goffset);
        }
        printf("GRAPH\n");
        for (auto e = edges.begin(); e != edges.end(); ++e) {
            for (size_t ve = 0; ve < e->second.size(); ++ve) {
                printf("%2d --> %2d\n", e->first, e->second[ve]);
            }
        }
    });
}


