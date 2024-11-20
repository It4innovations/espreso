
#include "dualgraph.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

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
        v->second.kernel.offset = offset;
        offset += v->second.kernel.size;
    }
    local = vertices;
}

void DualGraph::setFromDomains(const DecompositionFETI *decomposition, const std::vector<int> &lMap)
{
    for (size_t i = 0; i < lMap.size(); ) {
        int domains = lMap[i + 1];
        for (int d1 = 0; d1 < domains; ++d1) {
            for (int d2 = 0; d2 < domains; ++d2) {
                edges[d1].push_back(d2);
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
                sBuffer[n].push_back(v->second.kernel.offset);
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
            vertices[rBuffer[n][i]].kernel.offset = rBuffer[n][i + 1];
            vertices[rBuffer[n][i]].kernel.size = rBuffer[n][i + 2];
        }
    }

    std::map<int, int> coffset;
    for (size_t i = 0, offset = 0; i < lMap.size(); ) {
        int lsize = lMap[i];
        int domains = lMap[i + 1];
        int last = -1;
        for (int d1 = 0; d1 < domains; ++d1) {
            if (lMap[i + 2 + d1] < decomposition->dbegin) {
                int n = decomposition->noffset(lMap[i + 2 + d1]);
                if (last < n) {
                    for (int d2 = d1; d2 < domains; ++d2) {
                        if (decomposition->ismy(lMap[i + 2 + d2])) {
                            vertices[lMap[i + 2 + d2]].lambdas.indices.push_back({ coffset[d2], lsize });
                            vertices[lMap[i + 2 + d2]].lambdas.total += lsize;
                        }
                    }
                    last = n;
                }
            }
            if (decomposition->dend <= lMap[i + 2 + d1]) {
                vertices[lMap[i + 2 + d1]].lambdas.indices.push_back({ (int)offset, lsize });
                vertices[lMap[i + 2 + d1]].lambdas.total += lsize;
            }
            if (decomposition->ismy(lMap[i + 2 + d1])) {
                coffset[d1] += lsize;
            }
        }
        i += lMap[i + 1] + 2;
        offset += lsize;
    }

    for (auto v = vertices.begin(); v != vertices.end(); ++v) {
        printf("%2d :: %4d -> %4d\n", v->first, v->second.kernel.offset, v->second.kernel.size);
        for (size_t i = 0; i < v->second.lambdas.indices.size(); ++i) {
            printf("  %2d --> %2d\n", v->second.lambdas.indices[i].offset, v->second.lambdas.indices[i].size);
        }
    }
    for (auto e = edges.begin(); e != edges.end(); ++e) {
        for (size_t ve = 0; ve < edges[e->first].size(); ++ve) {
            printf("%2d --> %2d\n", e->first, edges[e->first][ve]);
        }
    }
}

void DualGraph::setFromClusters(const DecompositionFETI *decomposition, const std::vector<int> &lMap)
{

}
