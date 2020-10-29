
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
#include <iomanip>
#include <cfloat>
#include <fstream>

#include "basis/utilities/print.h"

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

void computeBodiesSurfacePlanes()
{
	if (info::mesh->surface->plane != NULL) {
		return;
	}
	profiler::syncstart("compute_surface_planes");

	std::vector<std::vector<Point> > planes(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		std::vector<Point> tplanes(2 * (info::mesh->surface->edistribution[t + 1] - info::mesh->surface->edistribution[t]));

		auto epointers = info::mesh->surface->epointers->datatarray().begin(t);
		auto enodes = info::mesh->surface->enodes->begin(t);
		for (size_t e = info::mesh->surface->edistribution[t], i = 0; e < info::mesh->surface->edistribution[t + 1]; ++e, ++i, ++epointers, ++enodes) {
			switch ((*epointers)->code) {
			case Element::CODE::LINE2:
			case Element::CODE::LINE3:
			{
				const Point &a = info::mesh->surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = info::mesh->surface->coordinates->datatarray()[enodes->at(1)];
				tplanes[2 * i] = (a + b) / 2;
				tplanes[2 * i + 1] = Point(b - a).normalize();
				std::swap(tplanes[2 * i + 1].x, tplanes[2 * i + 1].y);
				tplanes[2 * i + 1].y = -tplanes[2 * i + 1].y;
			} break;
			case Element::CODE::TRIANGLE3:
			case Element::CODE::TRIANGLE6:
			{
				const Point &a = info::mesh->surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = info::mesh->surface->coordinates->datatarray()[enodes->at(1)];
				const Point &c = info::mesh->surface->coordinates->datatarray()[enodes->at(2)];
				tplanes[2 * i] = (a + b + c) / 3;
				tplanes[2 * i + 1] = Point::cross(b - a, c - a).normalize();
			} break;
			case Element::CODE::SQUARE4:
			case Element::CODE::SQUARE8:
			{
				const Point &a = info::mesh->surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = info::mesh->surface->coordinates->datatarray()[enodes->at(1)];
				const Point &c = info::mesh->surface->coordinates->datatarray()[enodes->at(2)];
				const Point &d = info::mesh->surface->coordinates->datatarray()[enodes->at(3)];
				tplanes[2 * i] = (a + b + c + d) / 4;
				tplanes[2 * i + 1] = Point::cross(c - a, d - b).normalize();
			} break;
			default:
				eslog::error("ESPRESO internal error: unknown surface element.\n");
			}
		}

		planes[t].swap(tplanes);
	}

	info::mesh->surface->plane = new serializededata<esint, Point>(2, planes);

	DebugOutput::surfacePlanes("surface.planes", 1, 1);
	eslog::checkpointln("MESH: BODY SURFACE PLANES");
	profiler::syncend("compute_surface_planes");
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
			const auto &p = info::mesh->contacts->surface->plane->datatarray()[2 * esend[e]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&p), reinterpret_cast<const esint*>(&p) + 2 * sizeof(p) / sizeof(esint));
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

		// receive planes
		size = rBuffer[n][i];
		info::mesh->contacts->halo[n].plane = new serializededata<esint, Point>(2, tarray<Point>(info::env::OMP_NUM_THREADS, 2 * size));
		memcpy(info::mesh->contacts->halo[n].plane->datatarray().data(), rBuffer[n].data() + i + 1, 2 * sizeof(Point) * size);
		i += 1 + rBuffer[n][i] * (2 * sizeof(Point) / sizeof(esint));

		// receive global ids
		size = rBuffer[n][i];
		info::mesh->contacts->halo[n].nodes = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->halo[n].nodes->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(esint) * rBuffer[n][i]);
		i += 1 + rBuffer[n][i];

		// receive coordinates
		info::mesh->contacts->halo[n].coordinates = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->halo[n].coordinates->datatarray().data(), rBuffer[n].data() + i, sizeof(Point) * size);
	}

//	Communication::serialize([&] () {
//		std::cout << " -- " << info::mpi::rank << " --\n";
//		std::cout << info::mesh->contacts->surface->enodes->structures() << ": " << *info::mesh->contacts->surface->enodes << "\n";
//		std::cout << info::mesh->contacts->surface->coordinates->structures() << ": " << *info::mesh->contacts->surface->coordinates << "\n";
//
//		for (size_t n = 0; n < info::mesh->contacts->neighbors.size(); ++n) {
//			std::cout << info::mesh->contacts->halo[n].enodes->structures() << ": " << *info::mesh->contacts->halo[n].enodes << "\n";
//			std::cout << info::mesh->contacts->halo[n].coordinates->structures() << ": " << *info::mesh->contacts->halo[n].coordinates << "\n";
//		}
//	});

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
	std::vector<std::pair<esint, esint> > pair, atapair, nonself;

	auto o2n = [&] (esint offset) {
		return std::lower_bound(offsets.begin(), offsets.end(), offset + 1) - offsets.begin() - 1;
	};

	auto checkangle = [&] (esint i, esint j) {
		const Point *pi, *pj;
		if (i < offsets.front()) {
			pi = &(info::mesh->surface->plane->begin() + i)->at(1);
		} else {
			esint n = o2n(i);
			pi = &(info::mesh->contacts->halo[n].plane->begin() + (i - offsets[n]))->at(1);
		}
		if (j < offsets.front()) {
			pj = &(info::mesh->surface->plane->begin() + j)->at(1);
		} else {
			esint n = o2n(j);
			pj = &(info::mesh->contacts->halo[n].plane->begin() + (j - offsets[n]))->at(1);
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
//		if (
//				checkbody(i, j) &&
//				checkrange(i, j) &&
//				checkangle(i, j) &&
//				true) {
//
//			printf("%2dx%2d -> body: %d, range: %d[%f], angle: %d ", i, j, checkbody(i, j), checkrange(i, j), std::abs(eend[i].z - estart[j].z), checkangle(i, j));
//		}
		return
				checkbody(i, j) &&
				checkrange(i, j) &&
				checkangle(i, j) &&
				true;
	};

//	if (true) { // all-to-all
//		for (esint i = 0; i < offsets.front(); ++i) {
//			for (esint j = i + 1; j < offsets.back(); ++j) {
//				if (checkall(i, j)) {
//					if (i < j) {
//						atapair.push_back(std::make_pair(i, j));
//					} else {
//						atapair.push_back(std::make_pair(j, i));
//					}
//				}
//			}
//		}
////		printf("%d A-T-A AHIT: %lu\n", info::mpi::rank, atapair.size());
//	}

	if (true) { // tree
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
	} else {
		pair = atapair;
	}

//	printf("%3d [%d] HIT: %lu / %lu\n", info::mpi::rank, pair.size() == atapair.size(), pair.size(), atapair.size());

//	if (pair.size() != atapair.size()) {
//		printf("ERROR in global search\n");
//	}

	std::sort(pair.begin(), pair.end());
//	std::sort(atapair.begin(), atapair.end());
//	for (size_t i = 0, j = 0; i < atapair.size(); ++i) {
//		if (atapair[i] == pair[j]) {
//			++j;
//		} else {
//			printf("miss: %d -> %d\n", atapair[i].first, atapair[i].second);
//		}
//	}

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

//	Communication::serialize([&] () {
//		std::cout << " -- " << info::mpi::rank << " --\n";
//		std::cout << "local: " << *info::mesh->contacts->localPairs << "\n";
//		std::cout << "neigh: " << *info::mesh->contacts->neighPairs << "\n";
//	});

	DebugOutput::closeElements(1, 1);
	eslog::checkpointln("MESH: CLOSE ELEMENTS FOUND");
}

int intersection(const Point &p0, const Point &p1, const Point &q0, const Point &q1, Point &first, Point &second)
{
	double epsilon = 1e-16;
	auto perp = [] (Point &u, Point &v) {
		return u.x * v.y - u.y * v.x;
	};

	Point u = p1 - p0;
	Point v = q1 - q0;
	Point w = p0 - q0;
	double D = perp(u, v);

	// test if they are parallel (includes either being a point)
	if (fabs(D) < epsilon) { // lines are parallel
		if (perp(u, w) != 0 || perp(v, w) != 0) {
			return 0; // they are NOT collinear
		}
		// they are collinear

		// they are collinear segments - get  overlap (or not)
		double t0, t1; // endpoints of S1 in eqn for S2
		Point w2 = p1 - q0;
		if (v.x != 0) {
			t0 = w.x / v.x;
			t1 = w2.x / v.x;
		} else {
			t0 = w.y / v.y;
			t1 = w2.y / v.y;
		}
		if (t0 > t1) { // must have t0 smaller than t1
			double t = t0;
			t0 = t1;
			t1 = t;  // swap if not
		}
		if (t0 > 1 || t1 < 0) {
			return 0; // NO overlap
		}
		t0 = t0 < 0 ? 0 : t0; // clip to min 0
		t1 = t1 > 1 ? 1 : t1; // clip to max 1
		if (t0 == t1) { // intersect is a point
			first = q0 + v * t0;
			return 1;
		}

		// they overlap in a valid subsegment
		first = q0 + v * t0;
		first.z = t0;
		second = q0 + v * t1;
		second.z = t1;
		return 2;
	}

	// the segments are skew and may intersect in a point
	// get the intersect parameter for S1

	if (false) { // the second line is endless
		double s = perp(v, w) / D;
		if (s < 0 || s > 1) {
			return 0; // no intersect with S1
		}
	}

	// get the intersect parameter for S2
	double t = perp(u, w) / D;
	if (t < 0 || t > 1) {
		return 0; // no intersect with S2
	}

	first = q0 + v * t; // compute S1 intersect point
	first.z = t;
	return 1;
}

Point crossPoint(const Point &p0, const Point &p1, const Point &q0, const Point &q1)
{
	Point u = p1 - p0, v = q1 - q0, w = p0 - q0;
	double t = Point::cross2d(u, w) / Point::cross2d(u, v);
	return q0 + v * t;
}

int isIn(const Point &p, const std::vector<Point> &polygon)
{
	int n = polygon.size(), cn = 0;    // the crossing number counter

	for (size_t i = 0; i < polygon.size(); i++) {
		if (
				(polygon[i].y <= p.y && polygon[(i + 1) % n].y >  p.y) || // an upward crossing
				(polygon[i].y >  p.y && polygon[(i + 1) % n].y <= p.y)    // a downward crossing
			) {
			// compute the actual edge-ray intersect x-coordinate
			double vt = (p.y - polygon[i].y) / (polygon[(i + 1) % n].y - polygon[i].y);
			if (p.x < polygon[i].x + vt * (polygon[(i + 1) % n].x - polygon[i].x)) {
				++cn;   // a valid crossing of y = P.y right of P.x
			}
		}
	}
	return (cn & 1);    // 0 if even (out), and 1 if odd (in)
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

	esint print = 11;

	for (size_t i = 0; i < trias.size(); i += 3) {
		std::vector<Point> in(q), nextin;
		auto curr = [&] (size_t i) { return (i                ) % in.size(); };
		auto next = [&] (size_t i) { return (i + 1            ) % in.size(); };
		auto prev = [&] (size_t i) { return (i - 1 + in.size()) % in.size(); };

		Point p[3] = { trias[i], trias[i + 1], trias[i + 2] };
		Point line[3] = { (p[1] - p[0]).normalize(), (p[2] - p[1]).normalize(), (p[0] - p[2]).normalize() };

		if (id == print) printf(" :: [%5.2f,%5.2f][%5.2f,%5.2f][%5.2f,%5.2f] :: \n", p[0].x, p[0].y, p[1].x, p[1].y, p[2].x, p[2].y);
		for (int c = 0; c < 3; ++c) {
			std::vector<status> status(in.size());
			for (size_t j = 0; j < in.size(); ++j) {
				double distance = Point::cross2d(line[c], in[j] - p[c]);
				status[j] = eps < distance ? status::in : distance < -eps ? status::out : status::on;
				if (id == print) {
					printf("[%5.2f,%5.2f]", in[j].x, in[j].y);
					switch (status[j]) {
					case status::in: printf(": in\n"); break;
					case status::out: printf(": out\n"); break;
					case status::on: printf(": on\n"); break;
					case status::processed: printf(": processed\n"); break;
					}
				}
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
					if (id == print) printf("add cross: [%5.2f,%5.2f]\n", nextin.back().x, nextin.back().y);
					for (size_t k = 0; k < in.size(); ++k) {
						if (status[curr(j + k)] == status::in || status[curr(j + k)] == status::on) {
							nextin.push_back(in[curr(j + k)]); status[curr(j + k)] = status::processed;
							if (id == print) printf("add in[%lu]: [%5.2f,%5.2f]\n", curr(j + k), nextin.back().x, nextin.back().y);
						} else {
							nextin.push_back(crossPoint(p[c], p[(c + 1) % 3], in[prev(j + k)], in[curr(j + k)])); status[curr(j + k)] = status::processed;
							if (id == print) printf("add cross: [%5.2f,%5.2f]\n", nextin.back().x, nextin.back().y);
							break;
						}
					}
				}
			}
			if (id == print) {
				for (size_t j = 0; j < in.size(); ++j) {
					switch (status[j]) {
					case status::in: printf(": in\n"); break;
					case status::out: printf(": out\n"); break;
					case status::on: printf(": on\n"); break;
					case status::processed: printf(": processed\n"); break;
					}
				}
			}
			if (nextin.size()) {
				if (id == print) printf("crossed\n");
				in = nextin;
				nextin.clear();
			} else if (in.size() && status.front() != status::in) {
				if (id == print) printf("all-out\n");
				in.clear();
			}
			if (id == print) printf(" = = = = = = = = \n");
		}
		if (in.size()) {
			res.push_back(in);
		}
	}
	if (id == print) {
		for (size_t i = 0; i < res.size(); ++i) {
			printf("res: ");
			for (size_t j = 0; j < res[i].size(); ++j) {
				printf("[%5.2f,%5.2f]", res[i][j].x, res[i][j].y);
			}
			printf("\n");
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
	auto plane = info::mesh->surface->plane->cbegin();
	auto enodes = info::mesh->contacts->surface->enodes->cbegin();
	auto epointer = info::mesh->contacts->surface->epointers->datatarray();
	for (esint e = 0; e < info::mesh->surface->size; ++e, ++local, ++neigh, ++plane, ++enodes) {
		if (local->size() == 0 && neigh->size() == 0) {
			continue;
		}

		std::vector<Point> source, results;
		std::vector<std::vector<Point> > adepts;

		auto planee = info::mesh->surface->plane->begin() + e;
		const Point &pe = planee->at(0);
		const Point &ne = planee->at(1);
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
		source = ptrias;
		double fullArea = 0, piecedArea = 0;
		for (size_t i = 0; i < ptrias.size(); i += 3) {
			fullArea += .5 * Point::cross2d(ptrias[i + 1] - ptrias[i + 0], ptrias[i + 2] - ptrias[i + 0]);
		}

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
//				triangulate(p1, p2, intersection);
				adepts.push_back(p2);
				clip(e, ptrias, p2, intersection);
//				dummyClip(p1, p2, intersection);
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
					for (size_t j = 2; j < intersection[i].size(); j++) {
						results.push_back(intersection[i][0]);
						results.push_back(intersection[i][j - 1]);
						results.push_back(intersection[i][j]);
						piecedArea += .5 * Point::cross2d(intersection[i][j - 1] - intersection[i][0], intersection[i][j] - intersection[i][0]);
					}
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
//				triangulate(p1, p2, intersection);
//				dummyClip(p1, p2, intersection);
				adepts.push_back(p2);
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
					for (size_t j = 2; j < intersection[i].size(); j++) {
						results.push_back(intersection[i][0]);
						results.push_back(intersection[i][j - 1]);
						results.push_back(intersection[i][j]);
						piecedArea += .5 * Point::cross2d(intersection[i][j - 1] - intersection[i][0], intersection[i][j] - intersection[i][0]);
					}
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
		if (
				adepts.size()
				&& (fullArea <= 0.99 * piecedArea || 1.01 * piecedArea <= fullArea)
			) {
			printf("area[%d][%d]: %.10f != %.10f\n", info::mpi::rank, e, fullArea, piecedArea);
			std::ofstream os("clip" + std::to_string(info::mpi::rank) + "_" + std::to_string(e) + ".vtk");
			os << "# vtk DataFile Version 2.0\n";
			os << "EXAMPLE\n";
			os << "ASCII\n";
			os << "DATASET UNSTRUCTURED_GRID\n\n";
			esint points = source.size() + results.size();
			for (size_t aa = 0; aa < adepts.size(); ++aa) {
				points += adepts[aa].size();
			}
			os << "POINTS " << points << " float\n";
			for (size_t pp = 0; pp < source.size(); ++pp) {
				source[pp].z = 0;
				os << source[pp].x << " " << source[pp].y << " " << source[pp].z << "\n";
			}
			for (size_t aa = 0; aa < adepts.size(); ++aa) {
				for (size_t pp = 0; pp < adepts[aa].size(); ++pp) {
					adepts[aa][pp].z = 0;
					os << adepts[aa][pp].x << " " << adepts[aa][pp].y << " " << adepts[aa][pp].z << "\n";
				}
			}
			for (size_t pp = 0; pp < results.size(); ++pp) {
				results[pp].z = 0;
				os << results[pp].x << " " << results[pp].y << " " << results[pp].z << "\n";
			}

			esint cells = source.size() / 3 + results.size() / 3 + adepts.size();
			os << "\nCELLS " << cells << " " << cells + points << "\n";
			esint offset = 0;
			for (size_t pp = 0; pp < source.size(); pp += 3, offset += 3) {
				os << "3 " << offset << " " << offset + 1 << " " << offset + 2 << "\n";
			}
			for (size_t aa = 0; aa < adepts.size(); ++aa) {
				os << adepts[aa].size();
				for (size_t pp = 0; pp < adepts[aa].size(); ++pp) {
					os << " " << offset++;
				}
				os << "\n";
			}
			for (size_t pp = 0; pp < results.size(); pp += 3, offset += 3) {
				os << "3 " << offset << " " << offset + 1 << " " << offset + 2 << "\n";
			}

			os << "\nCELL_TYPES " << cells << "\n";
			for (size_t pp = 0; pp < source.size(); pp += 3) {
				os << "5\n";
			}
			for (size_t aa = 0; aa < adepts.size(); ++aa) {
				if (adepts[aa].size() == 3) {
					os << "5\n";
				}
				if (adepts[aa].size() == 4) {
					os << "9\n";
				}
			}
			for (size_t pp = 0; pp < results.size(); pp += 3) {
				os << "5\n";
			}
		}
	}

	info::mesh->contacts->interface = new serializededata<esint, esint>(dist, data);
	info::mesh->contacts->planeData = new serializededata<esint, double>(2, planedata);
	info::mesh->contacts->intersections = new serializededata<esint, Point>(3, triangles);

//	Communication::serialize([&] () {
//		std::cout << "interface: " << *info::mesh->contacts->interface << "\n";
//	});
//	esint iii = 0;
//	for (auto pp = info::mesh->contacts->planeData->begin(); pp != info::mesh->contacts->planeData->end(); ++pp) {
//		std::cout << iii++ << ": " << *pp << "\n";
//	}
//	std::cout << "plane: " << *info::mesh->contacts->planeData << "\n";

	DebugOutput::contact(1, 1);
	profiler::syncend("compute_contact_interface");
	eslog::checkpointln("MESH: CONTACT INTERFACE COMPUTED");
}

}
}
