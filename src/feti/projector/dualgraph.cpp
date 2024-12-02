
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

//    for (auto v = vertices.begin(); v != vertices.end(); ++v) {
////        printf("%2d :: %4d -> %4d\n", v->first, v->second.kernel.goffset, v->second.kernel.goffset + v->second.kernel.size);
//        printf("%2d, %2d ::", info::mpi::rank, v->first);
//        esint ll = 0;
//        for (size_t i = 0; i < v->second.lambdas.indices.size(); ++i) {
//            while (ll < v->second.lambdas.indices[i].offset) {
//                printf("    ");
//                ++ll;
//            }
//            for (esint j = 0; j < v->second.lambdas.indices[i].size; ++j, ++ll) {
//                printf(" %3d", v->second.lambdas.indices[i].offset + j);
//            }
//        }
//        printf("\n");
//    }
//    for (auto e = edges.begin(); e != edges.end(); ++e) {
//        for (size_t ve = 0; ve < edges[e->first].size(); ++ve) {
//            printf("%2d --> %2d\n", e->first, edges[e->first][ve]);
//        }
//    }
}

void DualGraph::setFromClusters(const DecompositionFETI *decomposition, const std::vector<int> &lMap)
{

}
