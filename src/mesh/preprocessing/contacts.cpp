
#include "meshpreprocessing.h"

#include "mesh/element.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/contactstore.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/structures/intervaltree.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/utils.h"
#include "basis/logging/timelogger.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "output/visualization/debug.h"

#include <algorithm>
#include <map>
#include <numeric>

namespace espreso {
namespace mesh {

void computeBodiesSurface()
{
	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), parents(threads), body(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointers(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = info::mesh->elements->procNodes->cbegin(t);
		auto neighs = info::mesh->elements->faceNeighbors->cbegin(t);
		const auto &epointers = info::mesh->elements->epointers->datatarray();

		std::vector<esint> fdist, fdata, fparents, fbody, ecounter((int)Element::CODE::SIZE);
		std::vector<Element*> fpointer;
		if (t == 0) {
			fdist.push_back(0);
		}

		for (size_t e = info::mesh->elements->distribution[t]; e < info::mesh->elements->distribution[t + 1]; ++e, ++neighs, ++nodes) {
			for (size_t n = 0; n < neighs->size(); ++n) {
				if (neighs->at(n) == -1) {
					auto face = epointers[e]->faces->begin() + n;
					for (auto f = face->begin(); f != face->end(); ++f) {
						fdata.push_back(nodes->at(*f));
					}
					fdist.push_back(fdata.size());
					fpointer.push_back(epointers[e]->facepointers->datatarray()[n]);
					fparents.push_back(e);
					fbody.push_back(info::mesh->elements->body->datatarray()[e]);
					++ecounter[(int)fpointer.back()->code];
				}
			}
		}

		facesDistribution[t].swap(fdist);

		faces[t].swap(fdata);
		fpointers[t].swap(fpointer);
		parents[t].swap(fparents);
		body[t].swap(fbody);
		ecounters[t].swap(ecounter);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t e = 0; e < ecounters[0].size(); e++) {
			ecounters[0][e] += ecounters[t][e];
		}
	}

	serializededata<esint, Element*>::balance(1, fpointers);
	info::mesh->surface->epointers = new serializededata<esint, Element*>(1, fpointers);
	info::mesh->surface->ecounters = ecounters[0];

	info::mesh->surface->edistribution = info::mesh->surface->epointers->datatarray().distribution();

	if (info::mesh->surface->ecounters[(int)Element::CODE::TRIANGLE3] == (esint)info::mesh->surface->edistribution.back()) {
		serializededata<esint, esint>::balance(3, faces, &info::mesh->surface->edistribution);
		info::mesh->surface->enodes = new serializededata<esint, esint>(3, faces);
		info::mesh->surface->triangles = info::mesh->surface->enodes;
		info::mesh->surface->tdistribution = info::mesh->surface->edistribution;
	} else {
		utils::threadDistributionToFullDistribution(facesDistribution);
		serializededata<esint, esint>::balance(facesDistribution, faces, &info::mesh->surface->edistribution);
		info::mesh->surface->enodes = new serializededata<esint, esint>(facesDistribution, faces);
	}
	serializededata<esint, esint>::balance(1, parents, &info::mesh->surface->edistribution);
	info::mesh->surface->parents = new serializededata<esint, esint>(1, parents);

	serializededata<esint, esint>::balance(1, body, &info::mesh->surface->edistribution);
	info::mesh->surface->body = new serializededata<esint, esint>(1, body);

	std::vector<esint> nodes(info::mesh->surface->enodes->datatarray().begin(), info::mesh->surface->enodes->datatarray().end());
	utils::sortAndRemoveDuplicates(nodes);
	info::mesh->surface->nodes = new serializededata<esint, esint>(1, tarray<esint>(threads, nodes));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = info::mesh->surface->enodes->datatarray().begin(t); n != info::mesh->surface->enodes->datatarray().end(t); ++n) {
			*n = std::lower_bound(info::mesh->surface->nodes->datatarray().begin(), info::mesh->surface->nodes->datatarray().end(), *n) - info::mesh->surface->nodes->datatarray().begin();
		}
	}

	std::vector<Point> coordinates;
	coordinates.reserve(nodes.size());
	for (size_t n = 0; n < nodes.size(); ++n) {
		coordinates.push_back(info::mesh->nodes->coordinates->datatarray()[nodes[n]]);
	}
	info::mesh->surface->coordinates = new serializededata<esint, Point>(1, tarray<Point>(info::mesh->surface->nodes->datatarray().distribution(), coordinates));

	info::mesh->surface->size = info::mesh->surface->edistribution.back();
	info::mesh->surface->offset = info::mesh->surface->edistribution.back();
	info::mesh->surface->totalSize = Communication::exscan(info::mesh->surface->offset);

	DebugOutput::surface("surface.bodies", 1, 1);
	eslog::checkpointln("MESH: BODY SURFACE COMPUTED");
}

void computeWarpedNormals(SurfaceStore * surface)
{
	if (surface->normal != NULL) {
		delete surface->normal;
	}
	if (surface->center != NULL) {
		delete surface->center;
	}
	profiler::syncstart("compute_warped_surface_normals");

	std::vector<std::vector<Point> > normal(info::env::OMP_NUM_THREADS), center(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		std::vector<Point> tnormal(surface->edistribution[t + 1] - info::mesh->surface->edistribution[t]);
		std::vector<Point> tcenter(surface->edistribution[t + 1] - info::mesh->surface->edistribution[t]);

		auto epointers = info::mesh->surface->epointers->datatarray().begin(t);
		auto enodes = surface->enodes->begin(t);
		for (size_t e = surface->edistribution[t], i = 0; e < surface->edistribution[t + 1]; ++e, ++i, ++epointers, ++enodes) {
			switch ((*epointers)->code) {
			case Element::CODE::LINE2:
			case Element::CODE::LINE3:
			{
				const Point &a = surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = surface->coordinates->datatarray()[enodes->at(1)];
				tcenter[i] = (a + b) / 2;
				tnormal[i] = Point(b - a).normalize();
				std::swap(tnormal[i].x, tnormal[i].y);
				tnormal[i].y = -tnormal[i].y;
			} break;
			case Element::CODE::TRIANGLE3:
			case Element::CODE::TRIANGLE6:
			{
				const Point &a = surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = surface->coordinates->datatarray()[enodes->at(1)];
				const Point &c = surface->coordinates->datatarray()[enodes->at(2)];
				tcenter[i] = (a + b + c) / 3;
				tnormal[i] = Point::cross(b - a, c - a).normalize();
			} break;
			case Element::CODE::SQUARE4:
			case Element::CODE::SQUARE8:
			{
				const Point &a = surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = surface->coordinates->datatarray()[enodes->at(1)];
				const Point &c = surface->coordinates->datatarray()[enodes->at(2)];
				const Point &d = surface->coordinates->datatarray()[enodes->at(3)];
				tcenter[i] = (a + b + c + d) / 4;
				tnormal[i] = Point::cross(c - a, d - b).normalize();
			} break;
			default:
				eslog::error("ESPRESO internal error: unknown surface element.\n");
			}
		}

		normal[t].swap(tnormal);
		center[t].swap(tcenter);
	}

	surface->normal = new serializededata<esint, Point>(1, normal);
	surface->center = new serializededata<esint, Point>(1, center);

	DebugOutput::warpedNormals("surface.planes", 1, 1);
	eslog::checkpointln("MESH: WARPED SURFACE NORMALS COMMPUTED");
	profiler::syncend("compute_warped_surface_normals");
}

void exchangeContactHalo()
{
	double area = info::ecf->input.contact.search_area;
	_Point<float> box[2];
	if (info::mesh->contacts->surface->coordinates->datatarray().size()) {
		box[0] = box[1] = info::mesh->nodes->coordinates->datatarray()[0];
	}
	for (auto c = info::mesh->contacts->surface->coordinates->datatarray().begin(); c < info::mesh->contacts->surface->coordinates->datatarray().end(); ++c) {
		box[0].x = std::min((float)c->x, box[0].x);
		box[0].y = std::min((float)c->y, box[0].y);
		box[0].z = std::min((float)c->z, box[0].z);
		box[1].x = std::max((float)c->x, box[1].x);
		box[1].y = std::max((float)c->y, box[1].y);
		box[1].z = std::max((float)c->z, box[1].z);
	}

	std::vector<_Point<float> > boxes(2 * info::mpi::size);
	Communication::allGather(box, boxes.data(), 6, MPI_FLOAT);

	auto areIntersected = [&] (_Point<float> *block1, _Point<float> *block2) {
		return	!
				(block1[1].x + area < block2[0].x || block2[1].x + area < block1[0].x) ||
				(block1[1].y + area < block2[0].y || block2[1].y + area < block1[0].y) ||
				(block1[1].z + area < block2[0].z || block2[1].z + area < block1[0].z);
	};

	info::mesh->contacts->neighbors.clear();
	for (int r = 0; r < info::mpi::size; r++) {
		if (r != info::mpi::rank) {
			if (areIntersected(box, boxes.data() + 2 * r)) {
				info::mesh->contacts->neighbors.push_back(r);
			}
		}
	}

	const auto &coordinates = info::mesh->contacts->surface->coordinates->datatarray();
	std::vector<std::vector<esint> > sBuffer(info::mesh->contacts->neighbors.size()), rBuffer(info::mesh->contacts->neighbors.size());
	for (size_t n = 0; n < info::mesh->contacts->neighbors.size(); ++n) {
		std::vector<esint> esend, nsend;
		const auto &min = boxes[2 * info::mesh->contacts->neighbors[n]];
		const auto &max = boxes[2 * info::mesh->contacts->neighbors[n] + 1];
		auto enodes = info::mesh->contacts->surface->enodes->begin();
		for (size_t e = 0; e < info::mesh->contacts->surface->enodes->structures(); ++e, ++enodes) {
			for (auto c = enodes->begin(); c != enodes->end(); ++c) {
				if (
						min.x - area <= coordinates[*c].x && coordinates[*c].x <= max.x + area &&
						min.y - area <= coordinates[*c].y && coordinates[*c].y <= max.y + area &&
						min.z - area <= coordinates[*c].z && coordinates[*c].z <= max.z + area) {

					esend.push_back(e);
					for (auto c = enodes->begin(); c != enodes->end(); ++c) {
						nsend.push_back(*c);
					}
					break;
				}
			}
		}
		size_t ssize = 0, nsize = nsend.size();
		ssize += 1 + 3 * esend.size(); // parent, body, epointer size
		ssize += 1 + esend.size() + nsend.size(); // enodes size
		ssize += 1 + esend.size() * (2 * sizeof(Point) / sizeof(esint)); // plane
		utils::sortAndRemoveDuplicates(nsend);
		ssize += 1 + nsend.size() * (1 + sizeof(Point) / sizeof(esint)); // ids + coordinates size

		sBuffer[n].reserve(ssize);

		// send parents and bodies
		sBuffer[n].push_back(esend.size());
		for (size_t e = 0; e < esend.size(); ++e) {
			sBuffer[n].push_back(info::mesh->contacts->surface->parents->datatarray()[esend[e]]);
		}
		for (size_t e = 0; e < esend.size(); ++e) {
			sBuffer[n].push_back(info::mesh->contacts->surface->body->datatarray()[esend[e]]);
		}
		for (size_t e = 0; e < esend.size(); ++e) {
			sBuffer[n].push_back(static_cast<int>(info::mesh->contacts->surface->epointers->datatarray()[esend[e]]->code));
		}

		// send enodes in target offsets
		sBuffer[n].push_back(esend.size() + nsize);
		enodes = info::mesh->contacts->surface->enodes->begin();
		for (size_t e = 0, prev = 0; e < esend.size(); prev = esend[e++]) {
			enodes += esend[e] - prev;
			sBuffer[n].push_back(enodes->size());
			for (auto c = enodes->begin(); c != enodes->end(); ++c) {
				sBuffer[n].push_back(std::lower_bound(nsend.begin(), nsend.end(), *c) - nsend.begin());
			}
		}

		// send planes
		sBuffer[n].push_back(esend.size());
		for (size_t e = 0; e < esend.size(); ++e) {
			const auto &nn = info::mesh->contacts->surface->normal->datatarray()[esend[e]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nn), reinterpret_cast<const esint*>(&nn) + sizeof(nn) / sizeof(esint));
		}
		for (size_t e = 0; e < esend.size(); ++e) {
			const auto &nc = info::mesh->contacts->surface->center->datatarray()[esend[e]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nc), reinterpret_cast<const esint*>(&nc) + sizeof(nc) / sizeof(esint));
		}

		// send global ids
		sBuffer[n].push_back(nsend.size());
		for (size_t c = 0; c < nsend.size(); ++c) {
			sBuffer[n].push_back(info::mesh->nodes->IDs->datatarray()[info::mesh->contacts->surface->nodes->datatarray()[nsend[c]]]);
		}

		// send coordinates
		for (size_t c = 0; c < nsend.size(); ++c) {
			const auto &p = info::mesh->contacts->surface->coordinates->datatarray()[nsend[c]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&p), reinterpret_cast<const esint*>(&p) + sizeof(p) / sizeof(esint));
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->contacts->neighbors)) {
		eslog::error("ESPRESO internal error: cannot exchange contact halo.\n");
	}

	info::mesh->contacts->halo.resize(info::mesh->contacts->neighbors.size());
	for (size_t n = 0, i = 0; n < info::mesh->contacts->neighbors.size(); ++n, i = 0) {
		// receive parents and bodies
		info::mesh->contacts->halo[n].parents = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
		memcpy(info::mesh->contacts->halo[n].parents->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(esint) * rBuffer[n][i]);
		info::mesh->contacts->halo[n].body = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
		memcpy(info::mesh->contacts->halo[n].body->datatarray().data(), rBuffer[n].data() + i + 1 + rBuffer[n][i], sizeof(esint) * rBuffer[n][i]);
		info::mesh->contacts->halo[n].epointers = new serializededata<esint, Element*>(1, tarray<Element*>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
		for (esint e = 0; e < rBuffer[n][i]; ++e) {
			info::mesh->contacts->halo[n].epointers->datatarray()[e] = &Mesh::edata[rBuffer[n][2 * rBuffer[n][i] + e + i + 1]];
		}
		i += 1 + 3 * rBuffer[n][i];

		// receive enodes
		std::vector<std::vector<esint> > edist(info::env::OMP_NUM_THREADS), enodes(info::env::OMP_NUM_THREADS);
		size_t size = i + rBuffer[n][i];
		++i;
		edist[0].push_back(0);
		while(i < size) {
			size_t nsize = rBuffer[n][i++];
			for (size_t nn = 0; nn < nsize; ++nn) {
				enodes[0].push_back(rBuffer[n][i++]);
			}
			edist[0].push_back(enodes[0].size());
		}
		serializededata<esint, esint>::balance(edist, enodes, &info::mesh->contacts->halo[n].parents->datatarray().distribution());
		info::mesh->contacts->halo[n].enodes = new serializededata<esint, esint>(edist, enodes);

		// receive warped normals
		size = rBuffer[n][i];
		info::mesh->contacts->halo[n].normal = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->halo[n].normal->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(Point) * size);
		i += 1 + size * (sizeof(Point) / sizeof(esint));

		info::mesh->contacts->halo[n].center = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->halo[n].center->datatarray().data(), rBuffer[n].data() + i, sizeof(Point) * size);
		i += size * (sizeof(Point) / sizeof(esint));

		// receive global ids
		size = rBuffer[n][i];
		info::mesh->contacts->halo[n].nodes = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->halo[n].nodes->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(esint) * rBuffer[n][i]);
		i += 1 + rBuffer[n][i];

		// receive coordinates
		info::mesh->contacts->halo[n].coordinates = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->halo[n].coordinates->datatarray().data(), rBuffer[n].data() + i, sizeof(Point) * size);
	}

	eslog::checkpointln("MESH: CLOSE BOUNDARY EXCHANGED");
}

void findCloseElements()
{
	std::vector<Point> estart, eend;

	auto mm = [&] (const serializededata<esint, esint> *enodes, const serializededata<esint, Point> *coo, size_t &offset) {
		for (auto e = enodes->begin(); e != enodes->end(); ++e, ++offset) {
			estart[offset] = eend[offset] = coo->datatarray()[e->front()];
			for (auto n = e->begin() + 1; n != e->end(); ++n) {
				estart[offset].x = std::min(estart[offset].x, coo->datatarray()[*n].x);
				estart[offset].y = std::min(estart[offset].y, coo->datatarray()[*n].y);
				estart[offset].z = std::min(estart[offset].z, coo->datatarray()[*n].z);
				eend[offset].x = std::max(eend[offset].x, coo->datatarray()[*n].x);
				eend[offset].y = std::max(eend[offset].y, coo->datatarray()[*n].y);
				eend[offset].z = std::max(eend[offset].z, coo->datatarray()[*n].z);
			}
		}
	};

	std::vector<esint> offsets;
	size_t offset = info::mesh->contacts->surface->enodes->structures();
	offsets.push_back(offset);
	for (size_t n = 0; n < info::mesh->contacts->neighbors.size(); ++n) {
		offset += info::mesh->contacts->halo[n].enodes->structures();
		offsets.push_back(offset);
	}
	estart.resize(offset);
	eend.resize(offset);
	offset = 0;
	mm(info::mesh->contacts->surface->enodes, info::mesh->contacts->surface->coordinates, offset);
	for (size_t n = 0; n < info::mesh->contacts->neighbors.size(); ++n) {
		mm(info::mesh->contacts->halo[n].enodes, info::mesh->contacts->halo[n].coordinates, offset);
	}

	IntervalTree tree(estart, eend);

	double eps = info::ecf->input.contact.search_area;
	double max_angle = -std::cos(M_PI * info::ecf->input.contact.max_angle / 180);
	bool self_contact = info::ecf->input.contact.self_contact;
	std::vector<std::pair<esint, esint> > pair;

	auto o2n = [&] (esint offset) {
		return std::lower_bound(offsets.begin(), offsets.end(), offset + 1) - offsets.begin() - 1;
	};

	auto checkangle = [&] (esint i, esint j) {
		const Point *pi, *pj;
		if (i < offsets.front()) {
			pi = &info::mesh->surface->normal->datatarray()[i];
		} else {
			esint n = o2n(i);
			pi = &info::mesh->contacts->halo[n].normal->datatarray()[i - offsets[n]];
		}
		if (j < offsets.front()) {
			pj = &info::mesh->surface->normal->datatarray()[j];
		} else {
			esint n = o2n(j);
			pj = &info::mesh->contacts->halo[n].normal->datatarray()[j - offsets[n]];
		}

		if (max_angle < (*pj) * (*pi)) { // the same direction < angle(PI / 2) == 0 < are opposite
			return false;
		}
		return true;
	};

	auto checkrange = [&] (esint i, esint j) {
		return !(
			eend[i].x + eps < estart[j].x || eend[j].x + eps < estart[i].x ||
			eend[i].y + eps < estart[j].y || eend[j].y + eps < estart[i].y ||
			eend[i].z + eps < estart[j].z || eend[j].z + eps < estart[i].z);
	};

	auto checkbody = [&] (esint i, esint j) {
		esint ibody, jbody;
		if (i < offsets.front()) {
			ibody = info::mesh->surface->body->datatarray()[i];
		} else {
			esint n = o2n(i);
			ibody = info::mesh->contacts->halo[n].body->datatarray()[i - offsets[n]];
		}
		if (j < offsets.front()) {
			jbody = info::mesh->surface->body->datatarray()[j];
		} else {
			esint n = o2n(j);
			jbody = info::mesh->contacts->halo[n].body->datatarray()[j - offsets[n]];
		}
		if (!self_contact && ibody == jbody) {
			return false;
		}
		if (i < offsets.front() && j < offsets.front()) {
			if (i < j) {
				return ibody < jbody;
			} else {
				return jbody < ibody;
			}
		}
		if (i < offsets.front()) {
			return ibody < jbody;
		}
		if (j < offsets.front()) {
			return jbody < ibody;
		}
		return false;
	};

	auto checkall = [&] (esint i, esint j) {
		return
				checkbody(i, j) &&
				checkrange(i, j) &&
				checkangle(i, j) &&
				true;
	};

	{ // traverse tree to find all close pairs
		pair.clear();
		for (esint i = std::exp2(tree.levels), first = i; i < std::exp2(tree.levels + 1); ++i) {
			esint begin = tree.begin(i);
			esint end = tree.end(i);
			if (begin == end) {
				continue;
			}

			Point min;
			tree.boxMin(i, min);

			auto check = [&] (esint p, esint begin, esint end) {
				for (auto pp = tree.permutation.cbegin() + begin; pp != tree.permutation.cbegin() + end; ++pp) {
					esint i = tree.permutation[p], j = *pp;
					if (checkall(i, j)) {
						if (i < j) {
							pair.push_back(std::make_pair(i, j));
						} else {
							pair.push_back(std::make_pair(j, i));
						}
					}
				}
			};

			std::function<void(size_t, size_t, esint)> traverse = [&] (size_t node, size_t max, esint p) {
				if (estart[tree.permutation[p]][tree.splitters[node].d] < tree.splitters[node].end + eps) {
					if (2 * node < tree.splitters.size()) {
						traverse(2 * node, max, p);
					} else {
						if (2 * node < max) {
							check(p, tree.begin(2 * node), tree.end(2 * node));
						}
					}
				}
				if (tree.splitters[node].start - eps < eend[tree.permutation[p]][tree.splitters[node].d]) {
					if (2 * node < tree.splitters.size()) {
						traverse(2 * node + 1, max, p);
					} else {
						if (2 * node + 1 < max) {
							check(p, tree.begin(2 * node + 1), tree.end(2 * node + 1));
						}
					}
				}
			};

			// go for all that are not in the same bucket
			if (tree.splitters.size() > 1) {
				for (auto p = tree.permutation.cbegin() + begin; p != tree.permutation.cbegin() + end; ++p) {
					if ((estart[*p].x <= min.x + eps) || (estart[*p].y <= min.y + eps) || (estart[*p].z <= min.z + eps)) {
						traverse(1, std::max(i, first), p - tree.permutation.cbegin());
					}
				}
			}

			// check faces within the bucket
			for (auto left = tree.permutation.cbegin() + begin; left != tree.permutation.cbegin() + end; ++left) {
				for (auto right = left + 1; right != tree.permutation.cbegin() + end; ++right) {
					esint i = *left, j = *right;
					if (checkall(i, j)) {
						if (i < j) {
							pair.push_back(std::make_pair(i, j));
						} else {
							pair.push_back(std::make_pair(j, i));
						}
					}
				}
			}
		}
	}

	std::sort(pair.begin(), pair.end());

	std::vector<esint> ldist = { 0 }, ldata, ndist = { 0 }, ndata;
	esint e = 0;
	for (size_t i = 0; i < pair.size() && pair[i].first < offsets.front(); ++i) {
		while (e < pair[i].first) {
			ldist.push_back(ldata.size());
			ndist.push_back(ndata.size());
			++e;
		}
		if (pair[i].second < offsets.front()) {
			ldata.push_back(pair[i].second);
		} else {
			ndata.push_back(o2n(pair[i].second));
			ndata.push_back(pair[i].second - offsets[ndata.back()]);
		}
	}
	while (e < offsets.front()) {
		ldist.push_back(ldata.size());
		ndist.push_back(ndata.size());
		++e;
	}

	info::mesh->contacts->localPairs = new serializededata<esint, esint>(ldist, ldata);
	info::mesh->contacts->neighPairs = new serializededata<esint, esint>(ndist, ndata);

	DebugOutput::closeElements(1, 1);
	eslog::checkpointln("MESH: CLOSE ELEMENTS FOUND");
}

void triangulate(std::vector<Point> &p, std::vector<Point> &triangles)
{
	auto prev= [&] (const int &i) -> int {
		return (i + p.size() - 1) % p.size();
	};
	auto next = [&] (const int &i) -> int {
		return (i + 1) % p.size();
	};

	if (p.size() == 3) {
		triangles = p;
		return;
	}
	if (p.size() == 4) {
		for (int v = 0, size = p.size(); v < size; ++v) {
			Point v1 = p[v] - p[prev(v)], v2 = p[next(v)] - p[prev(v)];
			if (Point::cross2d(v1.normalize(), v2.normalize()) < 0) {
				triangles.push_back(p[v]);
				triangles.push_back(p[next(v)]);
				triangles.push_back(p[next(next(v))]);
				triangles.push_back(p[next(next(v))]);
				triangles.push_back(p[prev(v)]);
				triangles.push_back(p[v]);
				return;
			}
		}
		if ((p[0] - p[2]).length() < (p[1] - p[3]).length()) {
			triangles.push_back(p[0]);
			triangles.push_back(p[1]);
			triangles.push_back(p[2]);
			triangles.push_back(p[0]);
			triangles.push_back(p[2]);
			triangles.push_back(p[3]);
		} else {
			triangles.push_back(p[0]);
			triangles.push_back(p[1]);
			triangles.push_back(p[3]);
			triangles.push_back(p[1]);
			triangles.push_back(p[2]);
			triangles.push_back(p[3]);
		}
		return;
	}
	eslog::error("ESPRESO internal error: cannot compute triangles from a polygon.\n");
}

void clip(esint id, const std::vector<Point> &trias, const std::vector<Point> &q, std::vector<std::vector<Point> > &res)
{
	double eps = 1e-10;
	enum status { in, out, on, processed };

	auto crossPoint = [] (const Point &p0, const Point &p1, const Point &q0, const Point &q1) {
		Point u = p1 - p0, v = q1 - q0, w = p0 - q0;
		double t = Point::cross2d(u, w) / Point::cross2d(u, v);
		return q0 + v * t;
	};

	for (size_t i = 0; i < trias.size(); i += 3) {
		std::vector<Point> in(q), nextin;
		auto curr = [&] (size_t i) { return (i                ) % in.size(); };
		auto next = [&] (size_t i) { return (i + 1            ) % in.size(); };
		auto prev = [&] (size_t i) { return (i - 1 + in.size()) % in.size(); };

		Point p[3] = { trias[i], trias[i + 1], trias[i + 2] };
		Point line[3] = { (p[1] - p[0]).normalize(), (p[2] - p[1]).normalize(), (p[0] - p[2]).normalize() };

		for (int c = 0; c < 3; ++c) {
			std::vector<status> status(in.size());
			for (size_t j = 0; j < in.size(); ++j) {
				double distance = Point::cross2d(line[c], in[j] - p[c]);
				status[j] = eps < distance ? status::in : distance < -eps ? status::out : status::on;
			}
			for (size_t j = 0; j < in.size(); ++j) {
				if (status[j] == status::in && status[next(j)] == status::on) {
					while (status[j = next(j)] == status::on) {
						status[j] = status::in;
					}
				}
			}
			for (size_t j = 0; j < in.size(); ++j) {
				if (status[j] == status::in && status[prev(j)] == status::on) {
					while (status[j = prev(j)] == status::on) {
						status[j] = status::in;
					}
				}
			}
			for (size_t j = 0; j < in.size(); ++j) {
				if (status[j] == status::in && status[prev(j)] == status::out) {
					nextin.push_back(crossPoint(p[c], p[(c + 1) % 3], in[prev(j)], in[j])); status[prev(j)] = status::processed;
					for (size_t k = 0; k < in.size(); ++k) {
						if (status[curr(j + k)] == status::in || status[curr(j + k)] == status::on) {
							nextin.push_back(in[curr(j + k)]); status[curr(j + k)] = status::processed;
						} else {
							nextin.push_back(crossPoint(p[c], p[(c + 1) % 3], in[prev(j + k)], in[curr(j + k)])); status[curr(j + k)] = status::processed;
							break;
						}
					}
				}
			}
			if (nextin.size()) {
				in = nextin;
				nextin.clear();
			} else if (in.size() && status.front() != status::in) {
				in.clear();
			}
		}
		if (in.size()) {
			res.push_back(in);
		}
	}
}

void computeContactInterface()
{
	profiler::syncstart("compute_contact_interface");

	double eps = info::ecf->input.contact.search_area;

	std::vector<Point> p1(20), p2(20), ptrias(60), ppp;
	std::vector<std::vector<Point> > intersection(20);
	Point z(0, 0, 1);

	std::vector<Point> triangles;
	std::vector<double> planedata;
	std::vector<esint> dist = { 0 }, data;

	auto local = info::mesh->contacts->localPairs->cbegin();
	auto neigh = info::mesh->contacts->neighPairs->cbegin();
	auto enodes = info::mesh->contacts->surface->enodes->cbegin();
	auto epointer = info::mesh->contacts->surface->epointers->datatarray();
	for (esint e = 0; e < info::mesh->surface->size; ++e, ++local, ++neigh, ++enodes) {
		if (local->size() == 0 && neigh->size() == 0) {
			continue;
		}

		const Point &pe = info::mesh->surface->center->datatarray()[e];
		const Point &ne = info::mesh->surface->normal->datatarray()[e];
		Point axis(0, 0, 1);
		double sin = 0, cos = 1;
		if (std::fabs(ne.z) < 0.999) {
			axis = Point::cross(ne, z).normalize();
			double angle = std::acos(ne.z);
			sin = std::sin(angle); cos = std::cos(angle);
		}
		double emin = pe.z, emax = pe.z;

		p1.clear();
		for (auto n = epointer[e]->polygon->begin(); n != epointer[e]->polygon->end(); ++n) {
			p1.push_back(info::mesh->contacts->surface->coordinates->datatarray()[enodes->at(*n)]);
			p1.back() -= pe;
			p1.back().rodrigues(axis, cos, sin);
			emin = std::min(emin, p1.back().z);
			emax = std::max(emax, p1.back().z);
		}
		ptrias.clear();
		triangulate(p1, ptrias);

		esint clips = 0;
		for (auto other = local->begin(); other != local->end(); ++other) {
			p2.clear();
			double min = 1e10, max = -1e10;
			auto onodes = info::mesh->contacts->surface->enodes->cbegin() + *other;
			for (auto n = epointer[*other]->polygon->rbegin(); n != epointer[*other]->polygon->rend(); ++n) {
				p2.push_back(info::mesh->contacts->surface->coordinates->datatarray()[onodes->at(*n)]);
				p2.back() -= pe;
				p2.back().rodrigues(axis, cos, sin);
				min = std::min(min, p2.back().z);
				max = std::max(max, p2.back().z);
			}
			if (!(max < emin - eps || emax + eps < min)) {
				intersection.clear();
				clip(e, ptrias, p2, intersection);
				if (intersection.size()) {
					if (clips == 0) {
						data.push_back(e);
						data.push_back(planedata.size());
						data.push_back(0);
						data.push_back(0);
						for (auto pp = p1.begin(); pp != p1.end(); ++pp) {
							planedata.push_back(pp->x);
							planedata.push_back(pp->y);
						}
					}
					++data[dist.back() + 2];
					data.push_back(*other);
					data.push_back(0);
					for (auto pp = p2.rbegin(); pp != p2.rend(); ++pp) {
						planedata.push_back(pp->x);
						planedata.push_back(pp->y);
					}
					for (size_t i = 0; i < intersection.size(); i++) {
						data.back() += intersection[i].size() - 2;
						for (size_t j = 2; j < intersection[i].size(); j++) {
							planedata.push_back(intersection[i][0].x);
							planedata.push_back(intersection[i][0].y);
							planedata.push_back(intersection[i][j - 1].x);
							planedata.push_back(intersection[i][j - 1].y);
							planedata.push_back(intersection[i][j].x);
							planedata.push_back(intersection[i][j].y);
						}
					}
					clips += intersection.size();
				}
				for (size_t i = 0; i < intersection.size(); i++) {
					for (size_t j = 0; j < intersection[i].size(); j++) {
						intersection[i][j].rodrigues(axis, cos, -sin);
						intersection[i][j] += pe;
					}
					for (size_t j = 2; j < intersection[i].size(); j++) {
						triangles.push_back(intersection[i][0]);
						triangles.push_back(intersection[i][j - 1]);
						triangles.push_back(intersection[i][j]);
					}
				}
			}
		}
		for (auto other = neigh->begin(); other != neigh->end(); ++other) {
			esint h = *other++;
			p2.clear();
			double min = 1e10, max = -1e10;
			auto onodes = info::mesh->contacts->halo[h].enodes->cbegin() + *other;
			auto oepointer = info::mesh->contacts->halo[h].epointers->datatarray();
			for (auto n = oepointer[*other]->polygon->rbegin(); n != oepointer[*other]->polygon->rend(); ++n) {
				p2.push_back(info::mesh->contacts->halo[h].coordinates->datatarray()[onodes->at(*n)]);
				p2.back() -= pe;
				p2.back().rodrigues(axis, cos, sin);
				min = std::min(min, p2.back().z);
				max = std::max(max, p2.back().z);
			}
			if (!(max < emin - eps || emax + eps < min)) {
				intersection.clear();
				clip(e, ptrias, p2, intersection);
				if (intersection.size()) {
					if (clips == 0) {
						data.push_back(e);
						data.push_back(planedata.size());
						data.push_back(0);
						data.push_back(0);
						for (auto pp = p1.begin(); pp != p1.end(); ++pp) {
							planedata.push_back(pp->x);
							planedata.push_back(pp->y);
						}
					}
					++data[dist.back() + 3];
					data.push_back(h);
					data.push_back(*other);
					data.push_back(0);
					for (auto pp = p2.rbegin(); pp != p2.rend(); ++pp) {
						planedata.push_back(pp->x);
						planedata.push_back(pp->y);
					}
					for (size_t i = 0; i < intersection.size(); i++) {
						data.back() += intersection[i].size() - 2;
						for (size_t j = 2; j < intersection[i].size(); j++) {
							planedata.push_back(intersection[i][0].x);
							planedata.push_back(intersection[i][0].y);
							planedata.push_back(intersection[i][j - 1].x);
							planedata.push_back(intersection[i][j - 1].y);
							planedata.push_back(intersection[i][j].x);
							planedata.push_back(intersection[i][j].y);
						}
					}
					clips += intersection.size();
				}
				for (size_t i = 0; i < intersection.size(); i++) {
					for (size_t j = 0; j < intersection[i].size(); j++) {
						intersection[i][j].rodrigues(axis, cos, -sin);
						intersection[i][j] += pe;
					}
					for (size_t j = 2; j < intersection[i].size(); j++) {
						triangles.push_back(intersection[i][0]);
						triangles.push_back(intersection[i][j - 1]);
						triangles.push_back(intersection[i][j]);
					}
				}
			}
		}
		if (clips) {
			dist.push_back(data.size());
		}
	}

	info::mesh->contacts->interface = new serializededata<esint, esint>(dist, data);
	info::mesh->contacts->planeData = new serializededata<esint, double>(2, planedata);
	info::mesh->contacts->intersections = new serializededata<esint, Point>(3, triangles);

	DebugOutput::contact(1, 1);
	profiler::syncend("compute_contact_interface");
	eslog::checkpointln("MESH: CONTACT INTERFACE COMPUTED");
}

}
}
