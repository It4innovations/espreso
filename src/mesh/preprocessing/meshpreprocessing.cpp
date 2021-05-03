
#include "meshpreprocessing.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/parser.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"

#include "mesh/element.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/bodystore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/surfacestore.h"
#include "math/matrix.dense.h"

#include "output/visualization/debug.h"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>

namespace espreso {
namespace mesh {

void computeBodies(NodeStore *nodes, ElementStore *elements, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors)
{
	profiler::syncstart("mesh_bodies_found");

	esint ebegin = elements->distribution.process.offset;
	esint eend = ebegin + elements->distribution.process.size;
	esint nbodies = 0, boffset;
	std::vector<int> body(elements->distribution.process.size, -1);

	{ // DFS for body search
		for (esint e = 0; e < elements->distribution.process.size; ++e) {
			std::vector<esint> stack;
			if (body[e] == -1) {
				stack.push_back(e);
				body[e] = nbodies;
				while (stack.size()) {
					esint current = stack.back();
					stack.pop_back();
					auto dual = elements->faceNeighbors->cbegin() + current;
					for (auto n = dual->begin(); n != dual->end(); ++n) {
						if (*n != -1 && ebegin <= *n && *n < eend) {
							if (body[*n - ebegin] == -1) {
								stack.push_back(*n - ebegin);
								body[*n - ebegin] = nbodies;
							}
						}
					}
				}
				nbodies++;
			}
		}

		// use unique indices across processes
		boffset = nbodies;
		Communication::exscan(boffset);
		for (size_t i = 0; i < body.size(); ++i) {
			body[i] += boffset;
		}
	}

	profiler::synccheckpoint("initial_dfs");

	std::unordered_map<esint, std::unordered_set<esint> > graph; // compact representation (vertex, edges)
	std::unordered_map<esint, esint> holders, labels;
	{ // exchange edges of the compact dual and compute compact representation
		for (esint b = 0; b < nbodies; ++b) {
			graph[b + boffset];
			holders[b + boffset] = info::mpi::rank;
		}
		struct ebody { esint e, b; };
		std::vector<esint> edistribution = Communication::getDistribution(elements->distribution.process.size);
		std::vector<std::vector<ebody> > sBuffer(neighbors.size()), rBuffer(neighbors.size());

		auto e2roffset = [&] (esint e) {
			auto rank = std::lower_bound(edistribution.begin(), edistribution.end(), e + 1) - edistribution.begin() - 1;
			return std::lower_bound(neighbors.begin(), neighbors.end(), rank) - neighbors.begin();
		};

		auto dual = elements->faceNeighbors->cbegin();
		for (esint e = 0; e < elements->distribution.process.size; ++e, ++dual) {
			for (auto n = dual->begin(); n != dual->end(); ++n) {
				if (*n != -1 && (*n < ebegin || eend <= *n)) {
					sBuffer[e2roffset(*n)].push_back(ebody{e + ebegin, body[e]});
				}
			}
		}

		for (size_t n = 0; n < sBuffer.size(); ++n) { // compute compact representation
			std::sort(sBuffer[n].begin(), sBuffer[n].end(), [] (const ebody &i, const ebody &j) { return i.e < j.e; });
			size_t i = 0, last = 0;
			while (i < sBuffer[n].size()) {
				last = i++;
				while (i < sBuffer[n].size() && sBuffer[n][last].b == sBuffer[n][i].b) {
					++i;
				}
				// i points to the first element with different body
				if (sBuffer[n][last].e != sBuffer[n][i - 1].e) {
					sBuffer[n][last] = sBuffer[n][i - 1];
				}
				++last;
			}
			sBuffer[n].resize(last);
		}

		if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighbors)) {
			eslog::internalFailure("cannot exchange bodies indices\n");
		}

		dual = elements->faceNeighbors->cbegin();
		for (esint e = 0; e < elements->distribution.process.size; ++e, ++dual) {
			for (auto n = dual->begin(); n != dual->end(); ++n) {
				if (*n != -1 && (*n < ebegin || eend <= *n)) {
					esint roffset = e2roffset(*n);
					esint nbody = std::lower_bound(rBuffer[roffset].begin(), rBuffer[roffset].end(), *n, [] (const ebody &ebody, esint eindex) { return ebody.e < eindex; })->b;
					graph[body[e]].insert(nbody);
					holders[nbody] = neighbors[roffset];
				}
			}
		}
	}

	profiler::synccheckpoint("graph_synchronization");

	std::unordered_set<int> nn;
	for (auto n = holders.begin(); n != holders.end(); ++n) {
		if (n->second != info::mpi::rank) {
			nn.insert(n->second);
		}
	}

	for (int groups = 1; groups < info::mpi::size; groups = groups << 1) {
		// remove vertices without edges (full body found)
		//////////////////////////////////////////////////
		for (auto v = graph.begin(); v != graph.end(); ) {
			if (v->second.empty()) {
				v = graph.erase(v);
			} else {
				++v;
			}
		}

		// compute the process that receives vertices without edge to other group
		/////////////////////////////////////////////////////////////////////////
		int sink = (info::mpi::rank / groups) % 2 == 0 ? info::mpi::rank + groups : info::mpi::rank - groups;
		if (info::mpi::size <= sink) {
			sink = info::mpi::rank - groups;
		}
		if (0 <= sink && sink < info::mpi::size) {
			nn.insert(sink);
			if (info::mpi::rank + groups < info::mpi::size && info::mpi::size <= info::mpi::rank + 2 * groups) {
				nn.insert(info::mpi::rank + groups);
			}
		} else {
			sink = info::mpi::rank;
		}

		std::vector<int> neighbors(nn.begin(), nn.end());
		std::sort(neighbors.begin(), neighbors.end());
		auto n2i = [&] (int n) {
			return std::lower_bound(neighbors.begin(), neighbors.end(), n) - neighbors.begin();
		};

		// compute target processes
		///////////////////////////
		std::unordered_map<esint, esint> targets = holders;
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				if ((info::mpi::rank / groups) % 2 != (targets[*n] / groups) % 2) {
					if (targets[v->first] == info::mpi::rank) {
						targets[v->first] = targets[*n];
					} else {
						if (std::abs(info::mpi::rank - targets[v->first]) < std::abs(info::mpi::rank - targets[*n])) {
							targets[v->first] = targets[*n];
						}
					}
				}
			}
		}

		// vertices that are without neighbor in the other half are sent to the sink process
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			if (targets[v->first] == info::mpi::rank) {
				targets[v->first] = sink;
			}
		}

		// serialize and exchange target processes
		//////////////////////////////////////////
		std::vector<esint> sTargets;
		std::vector<std::vector<esint> > rTargets(neighbors.size());
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			sTargets.push_back(v->first); sTargets.push_back(targets[v->first]);
		}
		if (!Communication::exchangeUnknownSize(sTargets, rTargets, neighbors)) {
			eslog::internalFailure("cannot exchange compact graph holders.\n");
		}
		for (size_t n = 0; n < rTargets.size(); ++n) {
			for (size_t i = 0; i < rTargets[n].size(); i += 2) {
				targets[rTargets[n][i]] = rTargets[n][i + 1];
			}
		}

		// serialize and exchange graph
		///////////////////////////////
		std::vector<std::vector<esint> > sGraph(neighbors.size()), rGraph(neighbors.size());
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			int t = n2i(targets[v->first]);
			sGraph[t].push_back(v->first);
			sGraph[t].push_back(v->second.size());
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				sGraph[t].push_back(*n);
				if ((targets[*n] / groups) % 2 == (neighbors[t] / groups) % 2) {
					sGraph[t].push_back(targets[*n]);
				} else {
					sGraph[t].push_back(holders[*n]);
				}
			}
		}

		// update holders of exchanged vertices
		///////////////////////////////////////
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				if ((targets[*n] / groups) % 2 == (info::mpi::rank / groups) % 2) {
					holders[*n] = targets[*n];
				}
			}
		}

		if (!Communication::exchangeUnknownSize(sGraph, rGraph, neighbors)) {
			eslog::internalFailure("cannot exchange compact graph.\n");
		}
		for (size_t n = 0; n < rGraph.size(); ++n) {
			for (size_t i = 0; i < rGraph[n].size(); ) {
				esint vertex = rGraph[n][i++];
				holders[vertex] = info::mpi::rank;
				for (esint e = 0; e < rGraph[n][i]; ++e) {
					graph[vertex].insert(rGraph[n][i + 1 + 2 * e]);
					holders[rGraph[n][i + 1 + 2 * e]] = rGraph[n][i + 1 + 2 * e + 1];
				}
				i += 2 * rGraph[n][i] + 1;
			}
		}

		// BFS to compute re-mapping (vertices are renumbered to the lowest number of a neighboring vertex)
		///////////////////////////////////////////////////////////////////////////////////////////////////
		std::unordered_map<esint, esint> remap;
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			esint min = v->first;
			std::unordered_set<esint> marked;
			if (remap.find(v->first) == remap.end()) {
				std::vector<esint> stack({ v->first });
				while (stack.size()) {
					auto e = graph[stack.back()];
					marked.insert(stack.back());
					stack.pop_back();
					for (auto n = e.begin(); n != e.end(); ++n) {
						if (graph.find(*n) != graph.end() && marked.find(*n) == marked.end()) {
							stack.push_back(*n);
							min = std::min(min, *n);
						}
					}
				}
			}
			for (auto m = marked.begin(); m != marked.end(); ++m) {
				remap[*m] = min;
				if (labels.find(*m) != labels.end() || (boffset <= *m && *m < boffset + nbodies)) {
					labels[*m] = min;
					labels[min] = min;
				}
			}
		}

		// exchange re-mapping data in order to synchronize edges
		/////////////////////////////////////////////////////////
		nn.clear();
		for (auto n = holders.begin(); n != holders.end(); ++n) {
			if (n->second != info::mpi::rank) {
				nn.insert(n->second);
			}
		}
		neighbors.assign(nn.begin(), nn.end());
		std::sort(neighbors.begin(), neighbors.end());

		std::vector<esint> sMap;
		std::vector<std::vector<esint> > rMap(neighbors.size());
		for (auto m = remap.begin(); m != remap.end(); ++m) {
			sMap.push_back(m->first); sMap.push_back(m->second);
		}
		if (!Communication::exchangeUnknownSize(sMap, rMap, neighbors)) {
			eslog::internalFailure("cannot exchange re-mapped vertices.\n");
		}
		for (size_t n = 0; n < rMap.size(); ++n) {
			for (size_t i = 0; i < rMap[n].size(); i += 2) {
				remap[rMap[n][i]] = rMap[n][i + 1];
			}
		}

		// update merged graph
		//////////////////////
		std::unordered_map<esint, std::unordered_set<esint> > mgraph;
		std::unordered_map<esint, esint> mholders;
		for (auto v = graph.begin(); v != graph.end(); ++v) {
			mholders[remap[v->first]] = info::mpi::rank;
			for (auto n = v->second.begin(); n != v->second.end(); ++n) {
				if (graph.find(*n) != graph.end()) {
					if (remap[v->first] != remap[*n]) {
						mgraph[remap[v->first]].insert(remap[*n]);
					}
				} else {
					mgraph[remap[v->first]].insert(remap[*n]);
					mholders[remap[*n]] = holders[*n];
				}
			}
		}

		mgraph.swap(graph);
		mholders.swap(holders);
	}

	for (esint b = boffset; b < boffset + nbodies; ++b) {
		auto it = labels.find(b);
		while (it != labels.end() && it->first != it->second) {
			it = labels.find(it->second);
		}
		labels[b] = it != labels.end() ? it->second : b;
	}

	bodies->size = 0;
	std::vector<esint> ulabels;
	for (esint b = boffset; b < boffset + nbodies; ++b) {
		if (labels[b] == b) {
			ulabels.push_back(b);
			++bodies->size;
		}
	}
	std::sort(ulabels.begin(), ulabels.end());
	bodies->offset = ulabels.size() ? ulabels.front() : 0;
	Communication::allGatherUnknownSize(ulabels);
	bodies->offset = std::lower_bound(ulabels.begin(), ulabels.end(), bodies->offset) - ulabels.begin();

	bodies->totalSize = ulabels.size();
	for (esint b = boffset; b < boffset + nbodies; ++b) {
		labels[b] = std::lower_bound(ulabels.begin(), ulabels.end(), labels[b]) - ulabels.begin();
	}
	for (size_t i = 0; i < body.size(); ++i) {
		body[i] = labels[body[i]];
	}

	if (elements->body == NULL) {
		elements->body = new serializededata<esint, int>(1, tarray<int>(elements->distribution.threads, body));
	} else {
		memcpy(elements->body->datatarray().data(), body.data(), elements->distribution.process.size * sizeof(int));
	}

	int rsize = bitMastSize(elementsRegions.size());
	std::vector<esint> bodyRegions(bodies->totalSize * rsize);
	for (esint e = 0; e < elements->distribution.process.size; ++e) {
		int b = elements->body->datatarray()[e];
		for (int r = 0; r < rsize; ++r) {
			bodyRegions[rsize * b + r] |= elements->regions->datatarray()[rsize * e + r];
		}
	}

	Communication::allReduce(bodyRegions.data(), NULL, bodyRegions.size(), MPITools::getType<esint>().mpitype, MPI_BOR);

	std::vector<esint> boffsets = { 0 };
	for (esint b = 0; b < bodies->totalSize; ++b) {
		int regions = 0;
		for (int r = 0; r < rsize; ++r) {
			for (size_t bit = 0; bit < 8 * sizeof(esint); ++bit) {
				if (bodyRegions[rsize * b + r] & ((esint)1 << bit)) {
					++regions;
				}
			}
		}
		boffsets.push_back(boffsets.back() + regions); // region ALL_ELEMENTS is not counted
	}

	for (int b = 0; b < bodies->totalSize; ++b) {
		for (int r = 0, rindex = 0; r < rsize; ++r) {
			for (size_t bit = 0; bit < 8 * sizeof(esint) && bit + r * sizeof(esint) < elementsRegions.size(); ++bit, ++rindex) {
				if (bodyRegions[b * rsize + r] & (1 << bit)) {
					elementsRegions[rindex]->bodies.push_back(b);
				}
			}
		}
	}

	std::vector<esint> bcount, esum(bodies->totalSize), fsum(bodies->totalSize);
	for (size_t r = 0; r < elementsRegions.size(); ++r) {
		std::fill(esum.begin(), esum.end(), 0);
		std::fill(fsum.begin(), fsum.end(), 0);
		for (auto e = elementsRegions[r]->elements->datatarray().begin(); e != elementsRegions[r]->elements->datatarray().end(); ++e) {
			int b = elements->body->datatarray()[*e];
			auto neighs = elements->faceNeighbors->begin() + *e;
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n == -1) {
					++fsum[b];
				}
			}
			++esum[b];
		}
		for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
			bcount.push_back(esum[elementsRegions[r]->bodies[b]]);
			bcount.push_back(fsum[elementsRegions[r]->bodies[b]]);
		}
	}

	Communication::allReduce(bcount, Communication::OP::SUM);

	for (size_t r = 0, offset = 0; r < elementsRegions.size(); ++r) {
		for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
			elementsRegions[r]->bodyElements.push_back(bcount[offset++]);
			elementsRegions[r]->bodyFaces.push_back(bcount[offset++]);
		}
	}

	profiler::syncend("mesh_bodies_found");
	eslog::checkpointln("MESH: MESH BODIES FOUND");
}

void linkNodesAndElements(ElementStore *elements, NodeStore *nodes, std::vector<int> &neighbors)
{
	linkNodesAndElements(
			nodes,
			neighbors,
			nodes->elements,
			elements->nodes,
			elements->IDs,
			elements->distribution.threads,
			true);
}

void linkNodesAndElements(
		NodeStore *nodes,
		std::vector<int> &neighbors,
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		std::vector<size_t> &edistribution,
		bool sortedIDs)
{
	profiler::syncstart("link_nodes_and_elements");
	size_t threads = info::env::OMP_NUM_THREADS;

	serializededata<esint, esint> *nIDs = nodes->IDs;
	serializededata<esint, int> *nranks = nodes->ranks;
	std::vector<size_t> &ndistribution = nodes->distribution;

	std::vector<esint> npermutation;
	if (!sortedIDs) {
		npermutation.resize(nIDs->datatarray().size());
		std::iota(npermutation.begin(), npermutation.end(), 0);
		std::sort(npermutation.begin(), npermutation.end(), [&] (esint i, esint j) {
			return nIDs->datatarray()[i] < nIDs->datatarray()[j];
		});
	}

	// thread x neighbor x vector(from, to)
	std::vector<std::vector<std::vector<std::pair<esint, esint> > > > sBuffer(threads);
	std::vector<std::vector<std::pair<esint, esint> > > rBuffer(neighbors.size());
	std::vector<std::pair<esint, esint> > localLinks;

	localLinks.resize(enodes->cend()->begin() - enodes->cbegin()->begin());
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto tnodes = enodes->cbegin(t);
		size_t offset = enodes->cbegin(t)->begin() - enodes->cbegin()->begin();

		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e, ++tnodes) {
			for (auto n = tnodes->begin(); n != tnodes->end(); ++n, ++offset) {
				localLinks[offset].first = *n;
				localLinks[offset].second = eIDs->datatarray()[e];
			}
		}
	}

	utils::sortWithInplaceMerge(localLinks, enodes->datatarray().distribution());
	profiler::synccheckpoint("local_links");

	std::vector<size_t> tbegin(threads);
	for (size_t t = 1; t < threads; t++) {
		tbegin[t] = std::lower_bound(localLinks.begin() + tbegin[t - 1], localLinks.end(), ndistribution[t], [] (std::pair<esint, esint> &p, esint n) { return p.first < n; }) - localLinks.begin();
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto ranks = nranks->cbegin(t);
		std::vector<std::vector<std::pair<esint, esint> > > tBuffer(neighbors.size());

		auto begin = localLinks.begin() + tbegin[t];
		auto end = begin;

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n, ++ranks) {
			while (begin != localLinks.end() && begin->first < (esint)n) ++begin;
			if (begin != localLinks.end() && begin->first == (esint)n) {
				end = begin;
				while (end != localLinks.end() && end->first == begin->first) ++end;
				esint nID = nIDs->datatarray()[begin->first];
				for (auto it = begin; it != end; ++it) {
					it->first = nID;
				}
				size_t i = 0;
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != info::mpi::rank) {
						while (neighbors[i] < *rank) ++i;
						tBuffer[i].insert(tBuffer[i].end(), begin, end);
					}
				}
				begin = end;
			}
		}

		sBuffer[t].swap(tBuffer);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t r = 0; r < sBuffer[0].size(); r++) {
			sBuffer[0][r].insert(sBuffer[0][r].end(), sBuffer[t][r].begin(), sBuffer[t][r].end());
		}
	}

	if (!sortedIDs) {
		std::sort(localLinks.begin(), localLinks.end());
		for (size_t n = 0; n < neighbors.size(); ++n) {
			std::sort(sBuffer[0][n].begin(), sBuffer[0][n].end());
		}
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, neighbors)) {
		eslog::internalFailure("addLinkFromTo - exchangeUnknownSize.\n");
	}
	profiler::synccheckpoint("exchange");

	std::vector<size_t> boundaries = { 0, localLinks.size() };
	for (size_t r = 0; r < rBuffer.size(); r++) {
		localLinks.insert(localLinks.end(), rBuffer[r].begin(), rBuffer[r].end());
		boundaries.push_back(localLinks.size());
	}

	utils::mergeAppendedData(localLinks, boundaries);
	if (!sortedIDs) {
		for (size_t i = 0, j = 0; i < localLinks.size(); i++) {
			while (nIDs->datatarray()[npermutation[j]] < localLinks[i].first) ++j;
			localLinks[i].first = npermutation[j];
		}
		std::sort(localLinks.begin(), localLinks.end());
	}

	std::vector<std::vector<esint> > linksBoundaries(threads);
	std::vector<std::vector<esint> > linksData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		if (ndistribution[t] != ndistribution[t + 1]) {
			auto llink = localLinks.begin();
			if (sortedIDs) {
				llink = std::lower_bound(localLinks.begin(), localLinks.end(), nIDs->datatarray()[ndistribution[t]], [] (std::pair<esint, esint> &p, esint n) { return p.first < n; });
			} else {
				llink = std::lower_bound(localLinks.begin(), localLinks.end(), ndistribution[t], [] (std::pair<esint, esint> &p, esint n) { return p.first < n; });
			}
			esint current;

			std::vector<esint> tBoundaries, tData;
			if (t == 0) {
				tBoundaries.push_back(0);
			}

			for (size_t n = ndistribution[t]; n < ndistribution[t + 1] && llink != localLinks.end(); ++n) {
				current = llink->first;
				if (
						(sortedIDs && current == nIDs->datatarray()[n]) ||
						(!sortedIDs && current == (esint)n)) {

					while (llink != localLinks.end() && current == llink->first) {
						tData.push_back(llink->second);
						++llink;
					}
				}
				tBoundaries.push_back(llink - localLinks.begin());
			}

			linksBoundaries[t].swap(tBoundaries);
			linksData[t].swap(tData);
		}
	}

	nelements = new serializededata<esint, esint>(linksBoundaries, linksData);

	profiler::synccheckpoint("rbuffer");
	profiler::syncend("link_nodes_and_elements");
	eslog::checkpointln("MESH: NODES AND ELEMENTS LINKED");
}

void computeNodesDuplication(NodeStore *nodes, std::vector<int> &neighborsWithMe)
{
	profiler::syncstart("compute_nodes_duplication");

	if (nodes->ranks) {
		delete nodes->ranks;
	}

	std::vector<esint> nids(nodes->IDs->datatarray().begin(), nodes->IDs->datatarray().end());

	std::vector<std::vector<esint> > sBuffer(neighborsWithMe.size(), nids), rBuffer(neighborsWithMe.size());
	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighborsWithMe)) {
		eslog::internalFailure("cannot exchange nodes ids.\n");
	}

	int threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > rdist(threads);
	std::vector<std::vector<int> >rdata(threads);
	rdist.front().push_back(0);
	std::vector<size_t> offset(neighborsWithMe.size());
	for (size_t n = 0; n < nids.size(); ++n) {
		for (size_t r = 0; r < neighborsWithMe.size(); ++r) {
			while (offset[r] < rBuffer[r].size() && rBuffer[r][offset[r]] < nids[n]) { ++offset[r]; }
			if (offset[r] < rBuffer[r].size() && rBuffer[r][offset[r]] == nids[n]) {
				rdata.front().push_back(neighborsWithMe[r]);
			}
		}
		rdist.front().push_back(rdata.front().size());
	}

	serializededata<esint, int>::balance(rdist, rdata);
	nodes->ranks = new serializededata<esint, int>(rdist, rdata);

	profiler::syncend("compute_nodes_duplication");
}

ElementStore* exchangeHalo(ElementStore *elements, NodeStore *nodes, std::vector<int> &neighbors)
{
	profiler::syncstart("exchange_halo");
	// halo elements are all elements that have some shared node

	esint ebegin = elements->distribution.process.offset;
	esint eend = ebegin + elements->distribution.process.size;

	size_t threads = info::env::OMP_NUM_THREADS;
	std::vector<std::vector<esint> > sBuffer(neighbors.size()), rBuffer(neighbors.size());

	std::vector<std::vector<std::vector<esint> > > hElements(threads);

	// we have to got through all nodes because intervals are not computed yet
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<esint> > telements(neighbors.size());
		auto elinks = nodes->elements->cbegin(t);
		size_t i = 0;

		for (auto ranks = nodes->ranks->cbegin(t); ranks != nodes->ranks->cend(t); ++ranks, ++elinks) {
			auto begin = elinks->begin();
			auto end = elinks->begin();
			if (ranks->size() > 1) {
				i = 0;
				while (begin != elinks->end() && *begin < ebegin) ++begin;
				end = begin;
				while (end != elinks->end() && *end < eend) ++end;
				for (auto rank = ranks->begin(); rank != ranks->end(); ++rank) {
					if (*rank != info::mpi::rank) {
						while (neighbors[i] < *rank) ++i;
						telements[i].insert(telements[i].end(), begin, end);
					}
				}
			}
		}
		hElements[t].swap(telements);
	}

	int rsize = bitMastSize(elements->regions->datatarray().size() / elements->distribution.process.size);

	std::vector<std::vector<size_t> > tdist(neighbors.size());
	for (size_t n = 0; n < neighbors.size(); ++n) {
		tdist[n] = { 0, hElements[0][n].size() };
	}
	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < neighbors.size(); ++n) {
			hElements[0][n].insert(hElements[0][n].end(), hElements[t][n].begin(), hElements[t][n].end());
			tdist[n].push_back(hElements[0][n].size());
		}
	}
	for (size_t n = 0; n < neighbors.size(); ++n) {
		utils::sortWithInplaceMerge(hElements[0][n], tdist[n]);
	}
	#pragma omp parallel for
	for (size_t n = 0; n < neighbors.size(); ++n) {
		utils::removeDuplicates(hElements[0][n]);
		tdist[n] = tarray<size_t>::distribute(threads, hElements[0][n].size());
		sBuffer[n].resize((4 + rsize) * hElements[0][n].size());
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const auto &IDs = elements->IDs->datatarray();
		const auto &body = elements->body->datatarray();
		const auto &material = elements->material->datatarray();
		const auto &epointer = elements->epointers->datatarray();
		const auto &regions = elements->regions->datatarray();
		for (size_t n = 0; n < neighbors.size(); ++n) {
			for (size_t e = tdist[n][t]; e < tdist[n][t + 1]; e++) {
				sBuffer[n][(4 + rsize) * e + 0] = IDs[hElements[0][n][e] - ebegin];
				sBuffer[n][(4 + rsize) * e + 1] = body[hElements[0][n][e] - ebegin];
				sBuffer[n][(4 + rsize) * e + 2] = material[hElements[0][n][e] - ebegin];
				sBuffer[n][(4 + rsize) * e + 3] = (esint)epointer[hElements[0][n][e] - ebegin]->code;
				memcpy(sBuffer[n].data() + (4 + rsize) * e + 4, regions.data() + (hElements[0][n][e] - ebegin) * rsize, sizeof(esint) * rsize);
			}
		}
	}
	profiler::synccheckpoint("sbuffer");

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighbors)) {
		eslog::internalFailure("exchange halo elements.\n");
	}
	profiler::synccheckpoint("exchange");

	std::vector<std::vector<esint> > hid(threads), hregions(threads);
	std::vector<std::vector<int> > hbody(threads), hmaterial(threads);

	std::vector<std::vector<Element*> > hcode(threads);

	for (size_t n = 0; n < rBuffer.size(); ++n) {
		std::vector<size_t> distribution = tarray<size_t>::distribute(threads, rBuffer[n].size() / (4 + rsize));
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; ++e) {
				hid[t].push_back(rBuffer[n][(4 + rsize) * e + 0]);
				hbody[t].push_back(rBuffer[n][(4 + rsize) * e + 1]);
				hmaterial[t].push_back(rBuffer[n][(4 + rsize) * e + 2]);
				hcode[t].push_back(&Mesh::edata[rBuffer[n][(4 + rsize) * e + 3]]);
				hregions[t].insert(hregions[t].end(), rBuffer[n].data() + (4 + rsize) * e + 4, rBuffer[n].data() + (4 + rsize) * e + 4 + rsize);
			}
		}
	}

	ElementStore *halo = new ElementStore();
	halo->IDs = new serializededata<esint, esint>(1, hid);
	halo->body = new serializededata<esint, int>(1, hbody);
	halo->material = new serializededata<esint, int>(1, hmaterial);
	halo->epointers = new serializededata<esint, Element*>(1, hcode);
	halo->regions = new serializededata<esint, esint>(rsize, hregions);

	halo->distribution.process.size = halo->IDs->datatarray().size();
	halo->distribution.threads = halo->IDs->datatarray().distribution();

	const auto &hIDs = halo->IDs->datatarray();
	std::vector<esint> permutation(hIDs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) { return hIDs[i] < hIDs[j]; });
	halo->permute(permutation);

	profiler::synccheckpoint("rbuffer");
	profiler::syncend("exchange_halo");
	eslog::checkpointln("MESH: HALO EXCHANGED");
	return halo;
}

void computeElementsFaceNeighbors(NodeStore *nodes, ElementStore *elements, std::vector<int> &neighbors)
{
	computeElementsNeighbors(
			nodes,
			neighbors,
			nodes->elements,
			elements->faceNeighbors,
			elements->nodes,
			elements->IDs,
			elements->epointers,
			elements->distribution.threads,
			[] (Element *e) { return e->faces; },
			false, // there are max 1 neighbor
			true); // sorted nodes IDs

	DebugOutput::faceNeighbors();
}

void computeElementsEdgeNeighbors(NodeStore *nodes, ElementStore *elements, std::vector<int> &neighbors)
{
	computeElementsNeighbors(
			nodes,
			neighbors,
			nodes->elements,
			elements->edgeNeighbors,
			elements->nodes,
			elements->IDs,
			elements->epointers,
			elements->distribution.threads,
			[] (Element *e) { return e->edges; },
			true, // we need to know the number of neighbors
			true); // sorted nodes IDs
}

void computeElementsNeighbors(
		NodeStore *nodes,
		std::vector<int> &neighbors,
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint>* &eneighbors,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		serializededata<esint, Element*> *epointers,
		std::vector<size_t> &edistribution,
		std::function<serializededata<int, int>*(Element*)> across,
		bool insertNeighSize,
		bool sortedIDs)
{
	profiler::syncstart("compute_elements_neighbors");
	if (nelements == NULL) {
		linkNodesAndElements(nodes, neighbors, nelements, enodes, eIDs, edistribution, sortedIDs);
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > dualDistribution(threads);
	std::vector<std::vector<esint> > dualData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = enodes->cbegin(t);

		std::vector<esint> tdist, tdata, intersection;
		if (t == 0) {
			tdist.push_back(0);
		}
		for (size_t e = edistribution[t]; e < edistribution[t + 1]; ++e, ++nodes) {
			bool hasNeighbor = false;
			for (auto interface = across(epointers->datatarray()[e])->begin(); interface != across(epointers->datatarray()[e])->end(); ++interface) {
				auto telements = nelements->cbegin() + nodes->at(*interface->begin());
				intersection.clear();
				for (auto n = telements->begin(); n != telements->end(); ++n) {
					if (*n != eIDs->datatarray()[e]) {
						intersection.push_back(*n);
					}
				}
				for (auto n = interface->begin() + 1; n != interface->end() && intersection.size(); ++n) {
					telements = nelements->cbegin() + nodes->at(*n);
					auto it1 = intersection.begin();
					auto it2 = telements->begin();
					auto last = intersection.begin();
					while (it1 != intersection.end()) {
						while (it2 != telements->end() && *it2 < *it1) {
							++it2;
						}
						if (it2 == telements->end()) {
							break;
						}
						if (*it1 == *it2) {
							*last++ = *it1++;
						} else {
							it1++;
						}
					}
					intersection.resize(last - intersection.begin());
				}
				if (insertNeighSize) {
					tdata.push_back(intersection.size());
					tdata.insert(tdata.end(), intersection.begin(), intersection.end());
				} else {
					tdata.push_back(intersection.size() ? intersection.front() : -1);
					if (intersection.size() > 1) {
						eslog::error("Input error: a face shared by 3 elements found.\n");
					}
				}
				hasNeighbor = hasNeighbor | intersection.size();
			}
			tdist.push_back(tdata.size());
//			if (!hasNeighbor && elements->distribution.size != 1) {
//				eslog::error("Input error: a dangling element found (the element without any neighbor).\n");
//			}
		}

		dualDistribution[t].swap(tdist);
		dualData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dualDistribution);

	eneighbors = new serializededata<esint, esint>(dualDistribution, dualData);

	profiler::syncend("compute_elements_neighbors");
	eslog::checkpointln("MESH: ELEMENTS NEIGHBOURS COMPUTED");
}

void computeSurfaceElementNeighbors(NodeStore *nodes, std::vector<int> &neigbors, SurfaceStore *surface)
{
	// out dated
	surface->IDs = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, surface->size));
	std::iota(surface->IDs->datatarray().begin(), surface->IDs->datatarray().end(), surface->offset);
	computeElementsNeighbors(
			nodes,
			neigbors,
			surface->nelements,
			surface->neighbors,
			surface->enodes,
			surface->IDs,
			surface->epointers,
			surface->edistribution,
			[] (Element *e) { return e->faces; },
			false, // there are max 1 neighbor
			false); // nodes IDs are not sorted
}

void computeElementsCenters(ElementStore *elements, NodeStore *nodes)
{
	if (elements->centers) {
		return;
	}
	profiler::syncstart("compute_element_centers");
	int threads = info::env::OMP_NUM_THREADS;

	elements->centers = new serializededata<esint, Point>(1, elements->distribution.threads);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		auto center = elements->centers->datatarray().begin(t);
		for (auto e = elements->nodes->cbegin(t); e != elements->nodes->cend(t); ++e, ++center) {
			for (auto n = e->begin(); n != e->end(); ++n) {
				*center += nodes->coordinates->datatarray()[*n];
			}
			*center /= e->size();
		}
	}
	profiler::syncend("compute_element_centers");
	eslog::checkpointln("MESH: ELEMENTS CENTERS COMPUTED");
}

void computeDecomposedDual(NodeStore *nodes, ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors, std::vector<esint> &dualDist, std::vector<esint> &dualData)
{
	profiler::syncstart("compute_decomposed_dual");
	bool separateRegions = info::ecf->input.decomposition.separate_regions;
	bool separateMaterials = info::ecf->input.decomposition.separate_materials;
	bool separateEtypes = info::ecf->input.decomposition.separate_etypes;

	size_t threads = info::env::OMP_NUM_THREADS;
	esint eBegin = elements->distribution.process.offset;
	esint eEnd   = eBegin + elements->distribution.process.size;

	std::vector<esint> dDistribution(elements->distribution.process.size + 1);
	std::vector<std::vector<esint> > dData(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tdata;
		int mat1 = 0, mat2 = 0, reg = 0, etype1 = 0, etype2 = 0;
		int rsize = bitMastSize(elementsRegions.size());

		auto neighs = elements->faceNeighbors->cbegin(t);
		for (size_t e = elements->distribution.threads[t]; e < elements->distribution.threads[t + 1]; ++e, ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1 && eBegin <= *n && *n < eEnd) {
					if (separateMaterials) {
						mat1 = elements->material->datatarray()[e];
						mat2 = elements->material->datatarray()[*n - eBegin];
					}
					if (separateRegions) {
						reg = memcmp(elements->regions->datatarray().data() + e * rsize, elements->regions->datatarray().data() + (*n - eBegin) * rsize, sizeof(esint) * rsize);
					}
					if (separateEtypes) {
						etype1 = (int)elements->epointers->datatarray()[e]->type;
						etype2 = (int)elements->epointers->datatarray()[*n - eBegin]->type;
					}

					if (mat1 == mat2 && !reg && etype1 == etype2) {
						tdata.push_back(*n - eBegin);
					}
				}
			}
			dDistribution[e + 1] = tdata.size();
		}

		dData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(dDistribution, elements->distribution.threads);
	for (size_t t = 1; t < threads; t++) {
		dData[0].insert(dData[0].end(), dData[t].begin(), dData[t].end());
	}

	dualDist.swap(dDistribution);
	dualData.swap(dData[0]);

	profiler::syncend("compute_decomposed_dual");
	eslog::checkpointln("MESH: LOCAL DUAL GRAPH COMPUTED");
}

void computeRegionsSurface(ElementStore *elements, NodeStore *nodes, ElementStore *halo, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors)
{
	profiler::syncstart("compute_region_surface");
	size_t threads = info::env::OMP_NUM_THREADS;
	esint eBegin = elements->distribution.process.offset;
	esint eEnd = eBegin + elements->distribution.process.size;

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
		std::vector<std::vector<Element*> > fpointers(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			esint hindex, addFace = 0;
			int rsize = bitMastSize(elementsRegions.size());
			auto nodes = elements->nodes->cbegin();
			auto neighs = elements->faceNeighbors->cbegin();
			const auto &regions = elements->regions->datatarray();
			const auto &epointers = elements->epointers->datatarray();

			std::vector<esint> fdist, fdata, ecounter((int)Element::CODE::SIZE);
			std::vector<Element*> fpointer;
			if (t == 0) {
				fdist.push_back(0);
			}

			esint prev = 0;
			for (auto e = elementsRegions[r]->elements->datatarray().cbegin(t); e != elementsRegions[r]->elements->datatarray().cend(t); prev = *e++) {
				nodes += *e - prev;
				neighs += *e - prev;
				for (size_t n = 0; n < neighs->size(); ++n) {
					if (neighs->at(n) != -1 && r) {
						if (neighs->at(n) < eBegin || eEnd <= neighs->at(n)) {
							hindex = std::lower_bound(halo->IDs->datatarray().begin(), halo->IDs->datatarray().end(), neighs->at(n)) - halo->IDs->datatarray().begin();
							addFace = memcmp(regions.data() + *e * rsize, halo->regions->datatarray().data() + hindex * rsize, sizeof(esint) * rsize);
						} else {
							addFace = memcmp(regions.data() + *e * rsize, regions.data() + (neighs->at(n) - eBegin) * rsize, sizeof(esint) * rsize);
						}
					} else {
						addFace = neighs->at(n) == -1;
					}
					if (addFace) {
						auto face = epointers[*e]->faces->begin() + n;
						for (auto f = face->begin(); f != face->end(); ++f) {
							fdata.push_back(nodes->at(*f));
						}
						fdist.push_back(fdata.size());
						fpointer.push_back(epointers[*e]->facepointers->datatarray()[n]);
						++ecounter[(int)fpointer.back()->code];
						addFace = 0;
					}
				}
			}

			facesDistribution[t].swap(fdist);
			faces[t].swap(fdata);
			fpointers[t].swap(fpointer);
			ecounters[t].swap(ecounter);
		}

		for (size_t t = 1; t < threads; t++) {
			for (size_t e = 0; e < ecounters[0].size(); e++) {
				ecounters[0][e] += ecounters[t][e];
			}
		}

		serializededata<esint, Element*>::balance(1, fpointers);
		elementsRegions[r]->surface->epointers = new serializededata<esint, Element*>(1, fpointers);
		elementsRegions[r]->surface->ecounters = ecounters[0];

		elementsRegions[r]->surface->edistribution = elementsRegions[r]->surface->epointers->datatarray().distribution();

		if (
				elementsRegions[r]->surface->edistribution.back() &&
				elementsRegions[r]->surface->ecounters[(int)Element::CODE::TRIANGLE3] == (esint)elementsRegions[r]->surface->edistribution.back()) {

			serializededata<esint, esint>::balance(3, faces, &elementsRegions[r]->surface->edistribution);
			elementsRegions[r]->surface->enodes = new serializededata<esint, esint>(3, faces);
			elementsRegions[r]->surface->triangles = elementsRegions[r]->surface->enodes;
			elementsRegions[r]->surface->tdistribution = elementsRegions[r]->surface->edistribution;
		} else {
			utils::threadDistributionToFullDistribution(facesDistribution);
			serializededata<esint, esint>::balance(facesDistribution, faces, &elementsRegions[r]->surface->edistribution);
			elementsRegions[r]->surface->enodes = new serializededata<esint, esint>(facesDistribution, faces);
		}
	}

	profiler::syncend("compute_region_surface");
	eslog::checkpointln("MESH: REGION SURFACE COMPUTED");
}

void triangularizeSurface(SurfaceStore *surface)
{
	if (surface == NULL) {
		return;
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	if (surface->triangles == NULL) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);


		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> ttriangles;
			std::vector<size_t> tintervals;
			if (t == 0) {
				tintervals.push_back(0);
			}

			auto elements = surface->enodes->cbegin(t);
			const auto &epointers = surface->epointers->datatarray().begin();

			for (size_t e = surface->edistribution[t]; e < surface->edistribution[t + 1]; ++e, ++elements) {
				for (auto n = epointers[e]->triangles->datatarray().cbegin(); n != epointers[e]->triangles->datatarray().cend(); ++n) {
					ttriangles.push_back(elements->at(*n));
				}
			}
			tintervals.push_back(ttriangles.size() / 3);

			intervals[t].swap(tintervals);
			triangles[t].swap(ttriangles);
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		surface->tdistribution = intervals[0];
		surface->triangles = new serializededata<esint, esint>(3, triangles);
	}

	eslog::checkpointln("MESH: SURFACE TRIANGULARIZED");
}

void triangularizeBoundary(BoundaryRegionStore *boundary)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	if (boundary->dimension == 2) {

		std::vector<std::vector<esint> > triangles(threads);
		std::vector<std::vector<size_t> > intervals(threads);


		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> ttriangles;
			std::vector<size_t> tintervals;
			if (t == 0) {
				tintervals.push_back(0);
			}

			auto elements = boundary->elements->cbegin(t);
			const auto &epointers = boundary->epointers->datatarray().begin();

			for (size_t e = boundary->distribution.threads[t]; e < boundary->distribution.threads[t + 1]; ++e, ++elements) {
				for (auto n = epointers[e]->triangles->datatarray().cbegin(); n != epointers[e]->triangles->datatarray().cend(); ++n) {
					ttriangles.push_back(elements->at(*n));
				}
			}
			tintervals.push_back(ttriangles.size() / 3);

			intervals[t].swap(tintervals);
			triangles[t].swap(ttriangles);
		}

		utils::threadDistributionToFullDistribution(intervals);
		utils::mergeThreadedUniqueData(intervals);

		boundary->triangles = new serializededata<esint, esint>(3, triangles);
	}

	eslog::checkpointln("MESH: BOUNDARY TRIANGULARIZED");
}

void computeBoundaryNodes(NodeStore *nodes, ElementStore *elements, DomainStore *domains, std::vector<int> &neighbors, std::vector<esint> &externalBoundary, std::vector<esint> &internalBoundary)
{
	profiler::syncstart("compute_boundary_nodes");
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > external(threads), internal(threads);

	esint eoffset = elements->distribution.process.offset;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> texternal, tinternal;

		auto neighbors = elements->faceNeighbors->cbegin() + domains->elements[domains->distribution[t]];
		auto enodes = elements->nodes->cbegin() + domains->elements[domains->distribution[t]];
		for (size_t d = domains->distribution[t]; d < domains->distribution[t + 1]; d++) {
			esint dbegin = domains->elements[d];
			esint dend = domains->elements[d + 1];

			for (esint e = dbegin; e < dend; ++e, ++neighbors, ++enodes) {
				auto epointer = elements->epointers->datatarray()[e];
				auto faces = epointer->faces->begin();

				for (size_t n = 0; n < neighbors->size(); ++n, ++faces) {
					if (neighbors->at(n) == -1) {
						for (auto f = faces->begin(); f != faces->end(); ++f) {
							texternal.push_back(enodes->at(*f));
						}
					} else if (neighbors->at(n) < dbegin + eoffset || dend + eoffset <= neighbors->at(n)) {
						for (auto f = faces->begin(); f != faces->end(); ++f) {
							tinternal.push_back(enodes->at(*f));
						}
					}
				}
			}
		}

		internal[t].swap(tinternal);
		external[t].swap(texternal);
	}

	utils::sortWithUniqueMerge(internal);
	utils::sortWithUniqueMerge(external);

	externalBoundary.swap(external[0]);

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(neighbors.begin(), neighbors.end(), neighbor) - neighbors.begin();
	};

	// external nodes need to be synchronized
	std::vector<std::vector<esint> > sBuffer(neighbors.size()), rBuffer(neighbors.size());
	std::vector<esint> nExternal;

	for (size_t i = 0; i < externalBoundary.size(); i++) {
		auto nrank = nodes->ranks->cbegin() + externalBoundary[i];
		for (auto rank = nrank->begin(); rank != nrank->end(); ++rank) {
			if (*rank != info::mpi::rank) {
				sBuffer[n2i(*rank)].push_back(nodes->IDs->datatarray()[externalBoundary[i]]);
			}
		}
	}

	for (size_t n = 0; n < neighbors.size(); n++) {
		std::sort(sBuffer[n].begin(), sBuffer[n].end());
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, neighbors)) {
		eslog::internalFailure("exchange external nodes.\n");
	}

	for (size_t n = 0; n < neighbors.size(); n++) {
		nExternal.insert(nExternal.end(), rBuffer[n].begin(), rBuffer[n].end());
	}
	utils::sortAndRemoveDuplicates(nExternal);

	for (size_t n = 0; n < neighbors.size(); n++) {
		nExternal.resize(std::set_difference(nExternal.begin(), nExternal.end(), sBuffer[n].begin(), sBuffer[n].end(), nExternal.begin()) - nExternal.begin());
	}

	for (size_t n = 0; n < nExternal.size(); n++) {
		std::vector<std::vector<esint> > tnExternal(threads);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			auto it = std::find(nodes->IDs->datatarray().cbegin() + nodes->distribution[t], nodes->IDs->datatarray().cbegin() + nodes->distribution[t + 1], nExternal[n]);
			if (it != nodes->IDs->datatarray().cbegin() + nodes->distribution[t + 1]) {
				tnExternal[t].push_back(it - nodes->IDs->datatarray().cbegin());
			}
		}

		for (size_t t = 0; t < threads; t++) {
			externalBoundary.insert(externalBoundary.end(), tnExternal[t].begin(), tnExternal[t].end());
		}
	}
	std::sort(externalBoundary.begin(), externalBoundary.end());

	internalBoundary.resize(internal[0].size());
	internalBoundary.resize(std::set_difference(internal[0].begin(), internal[0].end(), externalBoundary.begin(), externalBoundary.end(), internalBoundary.begin()) - internalBoundary.begin());

	profiler::syncend("compute_boundary_nodes");
	eslog::checkpointln("MESH: BOUNDARY NODES COMPUTED");
}

}
}
