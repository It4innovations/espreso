
#include "dualgraph.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

void DualGraph::clear()
{
    domains.vertices.clear();
    domains.edges.clear();
    clusters.vertices.clear();
    clusters.edges.clear();
}

void DualGraph::pushVertex(int offset, int index, int size)
{
    domains.vertices[index].offset = offset;
    domains.vertices[index].kernel.size = size;
}

void DualGraph::set(const DecompositionFETI *decomposition, const std::vector<int> &lMap)
{
    // init vertices
    int offset = 0;
    for (auto v = domains.vertices.begin(); v != domains.vertices.end(); ++v) {
        offset += v->second.kernel.size;
    }
    Communication::exscan(offset);
    for (auto v = domains.vertices.begin(); v != domains.vertices.end(); ++v) {
        v->second.rank = info::mpi::rank;
        v->second.kernel.goffset = offset;
        offset += v->second.kernel.size;
    }

    for (size_t i = 0; i < lMap.size(); ) {
        int ndomains = lMap[i + 1];
        for (int d1 = 0; d1 < ndomains; ++d1) {
            for (int d2 = 0; d2 < ndomains; ++d2) {
                domains.edges[lMap[i + 2 + d1]].push_back(lMap[i + 2 + d2]);
            }
        }
        i += lMap[i + 1] + 2;
    }
    for (auto e = domains.edges.begin(); e != domains.edges.end(); ++e) {
        utils::sortAndRemoveDuplicates(e->second);
    }

    { // synchronize domains
        std::vector<std::vector<int> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
        for (auto v = domains.vertices.begin(); v != domains.vertices.end(); ++v) {
            int last = -1;
            for (size_t ve = 0; ve < domains.edges[v->first].size(); ++ve) {
                int n = decomposition->noffset(domains.edges[v->first][ve]);
                if (!decomposition->ismy(domains.edges[v->first][ve]) && last < n) {
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
                domains.vertices[rBuffer[n][i]].rank = decomposition->neighbors[n];
                domains.vertices[rBuffer[n][i]].kernel.goffset = rBuffer[n][i + 1];
                domains.vertices[rBuffer[n][i]].kernel.size = rBuffer[n][i + 2];
            }
        }
    }

    for (size_t i = 0, offset = 0; i < lMap.size(); ) {
        int lsize = lMap[i];
        int ndomains = lMap[i + 1];
        for (int d = 0; d < ndomains; ++d) {
            domains.vertices[lMap[i + 2 + d]].lambdas.indices.push_back({ (int)offset, lsize });
            domains.vertices[lMap[i + 2 + d]].lambdas.total += lsize;
            std::vector<int> &lneighs = domains.vertices[lMap[i + 2 + d]].lambdas.indices.back().lower_neighs;
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

    // set clusters
    int csize = 0;
    for (auto v = domains.vertices.begin(); v != domains.vertices.end(); ++v) {
        csize = std::max(csize, v->second.kernel.size);
    }
    offset = csize;
    Communication::exscan(offset);
    for (auto v = domains.vertices.begin(); v != domains.vertices.end(); ++v) {
        clusters.vertices[v->second.rank].rank = v->second.rank;
        clusters.vertices[v->second.rank].offset.push_back(v->second.offset);
        clusters.vertices[v->second.rank].kernel.goffset = offset;
        clusters.vertices[v->second.rank].kernel.size = csize;
    }

    for (auto e = domains.edges.begin(); e != domains.edges.end(); ++e) {
        for (size_t i = 0; i < e->second.size(); ++i) {
            clusters.edges[domains.vertices[e->first].rank].push_back(domains.vertices[e->second[i]].rank);
        }
        utils::sortAndRemoveDuplicates(clusters.edges[domains.vertices[e->first].rank]);
    }

    { // synchronize clusters
        std::vector<std::vector<int> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
        for (auto v = clusters.vertices.begin(); v != clusters.vertices.end(); ++v) {
            if (v->second.rank != info::mpi::rank) continue;

            for (size_t ve = 0; ve < clusters.edges[v->first].size(); ++ve) {
                if (clusters.edges[v->first][ve] != info::mpi::rank) {
                    int n = std::lower_bound(decomposition->neighbors.begin(), decomposition->neighbors.end(), clusters.edges[v->first][ve]) - decomposition->neighbors.begin();
                    sBuffer[n].push_back(v->first);
                    sBuffer[n].push_back(v->second.kernel.goffset);
                    sBuffer[n].push_back(v->second.kernel.size);
                }
            }
        }

        if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, decomposition->neighbors)) {
            eslog::error("cannot exchange dual graph info\n");
        }

        for (size_t n = 0; n < rBuffer.size(); ++n) {
            for (size_t i = 0; i < rBuffer[n].size(); i += 3) {
                clusters.vertices[rBuffer[n][i]].rank = decomposition->neighbors[n];
                clusters.vertices[rBuffer[n][i]].kernel.goffset = rBuffer[n][i + 1];
                clusters.vertices[rBuffer[n][i]].kernel.size = rBuffer[n][i + 2];
            }
        }
    }

    for (size_t i = 0, offset = 0; i < lMap.size(); ) {
        int lsize = lMap[i];
        int ndomains = lMap[i + 1];
        for (int d = 0, prev = -1; d < ndomains; ++d) {
            int c = domains.vertices[lMap[i + 2 + d]].rank;
            if (prev == c) { continue; }
            prev = c;
            clusters.vertices[c].lambdas.indices.push_back({ (int)offset, lsize });
            clusters.vertices[c].lambdas.total += lsize;
            std::vector<int> &lneighs = clusters.vertices[c].lambdas.indices.back().lower_neighs;
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

void DualGraph::spread(const DecompositionFETI *decomposition)
{
    { // domains
        std::map<int, std::vector<int> > next = domains.edges;

        for (auto e = domains.edges.begin(); e != domains.edges.end(); ++e) {
            for (size_t ve = 0; ve < e->second.size(); ++ve) {
                for (size_t vt = 0; vt < domains.edges[e->second[ve]].size(); ++vt) {
                    next[e->first].push_back(domains.edges[e->second[ve]][vt]);
                }
            }
        }

        for (auto e = next.begin(); e != next.end(); ++e) {
            utils::sortAndRemoveDuplicates(e->second);
        }

        std::vector<std::vector<int> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
        if (decomposition->neighbors.size()) {
            sBuffer[0].push_back(domains.vertices.size());
            for (auto v = domains.vertices.begin(); v != domains.vertices.end(); ++v) {
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
                domains.vertices[rBuffer[n][i]].rank = rBuffer[n][i + 1];
                domains.vertices[rBuffer[n][i]].kernel.goffset = rBuffer[n][i + 2];
                domains.vertices[rBuffer[n][i]].kernel.size    = rBuffer[n][i + 3];
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
        domains.edges.swap(next);
    }

    { // clusters
        std::map<int, std::vector<int> > next = clusters.edges;

        for (auto e = clusters.edges.begin(); e != clusters.edges.end(); ++e) {
            for (size_t ve = 0; ve < e->second.size(); ++ve) {
                for (size_t vt = 0; vt < clusters.edges[e->second[ve]].size(); ++vt) {
                    next[e->first].push_back(clusters.edges[e->second[ve]][vt]);
                }
            }
        }

        for (auto e = next.begin(); e != next.end(); ++e) {
            utils::sortAndRemoveDuplicates(e->second);
        }

        std::vector<std::vector<int> > sBuffer(decomposition->neighbors.size()), rBuffer(decomposition->neighbors.size());
        if (decomposition->neighbors.size()) {
            sBuffer[0].push_back(clusters.vertices.size());
            for (auto v = clusters.vertices.begin(); v != clusters.vertices.end(); ++v) {
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
            eslog::error("cannot spread cluster dual graph.\n");
        }

        for (size_t n = 0; n < rBuffer.size(); ++n) {
            size_t i = 1;
            for (int v = 0; v < rBuffer[n][0]; ++v, i += 4) {
                clusters.vertices[rBuffer[n][i]].rank = rBuffer[n][i + 1];
                clusters.vertices[rBuffer[n][i]].kernel.goffset = rBuffer[n][i + 2];
                clusters.vertices[rBuffer[n][i]].kernel.size    = rBuffer[n][i + 3];
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
        clusters.edges.swap(next);
    }
}

void DualGraph::print()
{
    Communication::serialize([&] () {
        if (info::mpi::rank == 0) {
            printf(" >> DOMAINS << \n");
        }
        printf(" --- rank  %2d ---\n", info::mpi::rank);
        for (auto v = domains.vertices.begin(); v != domains.vertices.end(); ++v) {
            printf("%2d :: holder:%2d, goffset:%2d\n", v->first, v->second.rank, v->second.kernel.goffset);
        }
        printf("GRAPH\n");
        for (auto e = domains.edges.begin(); e != domains.edges.end(); ++e) {
            for (size_t ve = 0; ve < e->second.size(); ++ve) {
                printf("%2d --> %2d\n", e->first, e->second[ve]);
            }
        }
    });

    Communication::serialize([&] () {
        if (info::mpi::rank == 0) {
            printf(" >> CLUSTERS << \n");
        }
        printf(" --- rank  %2d ---\n", info::mpi::rank);
        for (auto v = clusters.vertices.begin(); v != clusters.vertices.end(); ++v) {
            printf("%2d :: holder:%2d, goffset:%2d\n", v->first, v->second.rank, v->second.kernel.goffset);
        }
        printf("GRAPH\n");
        for (auto e = clusters.edges.begin(); e != clusters.edges.end(); ++e) {
            for (size_t ve = 0; ve < e->second.size(); ++ve) {
                printf("%2d --> %2d\n", e->first, e->second[ve]);
            }
        }
    });
}


