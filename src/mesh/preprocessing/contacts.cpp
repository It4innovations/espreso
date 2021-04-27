
#include "meshpreprocessing.h"

#include "mesh/element.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/bodystore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/contactstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/contactinterfacestore.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/structures/intervaltree.h"
#include "wrappers/mpi/communication.h"
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
#include <unordered_set>
#include <unordered_map>

namespace espreso {
namespace mesh {

void computeBodiesSurface()
{
	for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
		switch (it->second.detection) {
		case ContactInterfaceConfiguration::DETECTION::ALL_BODIES:
			for (size_t r = 0; r < info::mesh->elementsRegions.size(); ++r) {
				info::mesh->elementsRegions[r]->contact.gap = std::max(it->second.gap, info::mesh->elementsRegions[r]->contact.gap);
				info::mesh->elementsRegions[r]->contact.angle = std::max(it->second.angle, info::mesh->elementsRegions[r]->contact.angle);
				info::mesh->elementsRegions[r]->contact.self_contact |= it->second.self_contact;
			}
			break;
		case ContactInterfaceConfiguration::DETECTION::BODY_LIST:
			for (size_t b = 0; b < it->second.body_list.size(); ++b) {
				info::mesh->eregion(it->second.body_list[b])->contact.gap = std::max(it->second.gap, info::mesh->eregion(it->second.body_list[b])->contact.gap);
				info::mesh->eregion(it->second.body_list[b])->contact.angle = std::max(it->second.angle, info::mesh->eregion(it->second.body_list[b])->contact.angle);
				info::mesh->eregion(it->second.body_list[b])->contact.self_contact |= it->second.self_contact;
			}
			break;
		case ContactInterfaceConfiguration::DETECTION::CONTACT_PAIR:
			eslog::internalFailure("implement CONTACT_PAIR detection.\n");
			break;
		}
	}

	info::mesh->elements->contact = new serializededata<esint, ContactInfo>(1, tarray<ContactInfo>(info::mesh->elements->epointers->datatarray().distribution(), 1, ContactInfo()));
	for (size_t r = 0; r < info::mesh->elementsRegions.size(); ++r) {
		for (auto e = info::mesh->elementsRegions[r]->elements->datatarray().begin(); e != info::mesh->elementsRegions[r]->elements->datatarray().end(); ++e) {
			info::mesh->elements->contact->datatarray()[*e].gap = std::max(info::mesh->elementsRegions[r]->contact.gap, info::mesh->elements->contact->datatarray()[*e].gap);
			info::mesh->elements->contact->datatarray()[*e].angle = std::max(info::mesh->elementsRegions[r]->contact.angle, info::mesh->elements->contact->datatarray()[*e].angle);
			info::mesh->elements->contact->datatarray()[*e].self_contact |= info::mesh->elementsRegions[r]->contact.self_contact;
		}
	}

	if (info::mesh->elements->faceNeighbors == NULL) {
		computeElementsFaceNeighbors();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), parents(threads), body(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
	std::vector<std::vector<Element*> > fpointers(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nodes = info::mesh->elements->nodes->cbegin(t);
		auto neighs = info::mesh->elements->faceNeighbors->cbegin(t);
		const auto &epointers = info::mesh->elements->epointers->datatarray();

		std::vector<esint> fdist, fdata, fparents, fbody, ecounter((int)Element::CODE::SIZE);
		std::vector<Element*> fpointer;
		if (t == 0) {
			fdist.push_back(0);
		}

		for (size_t e = info::mesh->elements->distribution.threads[t]; e < info::mesh->elements->distribution.threads[t + 1]; ++e, ++neighs, ++nodes) {
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
	info::mesh->surface->nIDs = new serializededata<esint, esint>(1, tarray<esint>(threads, nodes));

	std::vector<std::vector<ContactInfo> > contact(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = info::mesh->surface->enodes->datatarray().begin(t); n != info::mesh->surface->enodes->datatarray().end(t); ++n) {
			*n = std::lower_bound(info::mesh->surface->nodes->datatarray().begin(), info::mesh->surface->nodes->datatarray().end(), *n) - info::mesh->surface->nodes->datatarray().begin();
		}
		for (auto n = info::mesh->surface->nIDs->datatarray().begin(t); n != info::mesh->surface->nIDs->datatarray().end(t); ++n) {
			*n = info::mesh->nodes->IDs->datatarray()[*n];
		}
		std::vector<ContactInfo> tcontact;
		tcontact.reserve(info::mesh->surface->parents->datatarray().size(t));
		for (auto e = info::mesh->surface->parents->datatarray().begin(t); e != info::mesh->surface->parents->datatarray().end(t); ++e) {
			tcontact.push_back(info::mesh->elements->contact->datatarray()[*e]);
		}
		contact[t].swap(tcontact);
	}
	info::mesh->surface->contact = new serializededata<esint, ContactInfo>(1, contact);

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
	if (surface->base != NULL) {
		delete surface->base;
	}
	profiler::syncstart("compute_warped_surface_normals");

	std::vector<std::vector<Point> > normal(info::env::OMP_NUM_THREADS), base(info::env::OMP_NUM_THREADS), parameters(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		std::vector<Point> tnormal(surface->edistribution[t + 1] - info::mesh->surface->edistribution[t]);
		std::vector<Point> tparamters(2 * (surface->edistribution[t + 1] - info::mesh->surface->edistribution[t]));
		std::vector<Point> tbase(surface->edistribution[t + 1] - info::mesh->surface->edistribution[t]);

		auto epointers = info::mesh->surface->epointers->datatarray().begin(t);
		auto enodes = surface->enodes->begin(t);
		for (size_t e = surface->edistribution[t], i = 0; e < surface->edistribution[t + 1]; ++e, ++i, ++epointers, ++enodes) {
			switch ((*epointers)->code) {
//			case Element::CODE::LINE2:
//			case Element::CODE::LINE3:
//			{
//				const Point &a = surface->coordinates->datatarray()[enodes->at(0)];
//				const Point &b = surface->coordinates->datatarray()[enodes->at(1)];
//				tbase[i] = (a + b) / 2;
//				tnormal[i] = Point(b - a).normalize();
//				std::swap(tnormal[i].x, tnormal[i].y);
//				tnormal[i].y = -tnormal[i].y;
//			} break;
			case Element::CODE::TRIANGLE3:
//			case Element::CODE::TRIANGLE6:
			{
				const Point &a = surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = surface->coordinates->datatarray()[enodes->at(1)];
				const Point &c = surface->coordinates->datatarray()[enodes->at(2)];
				tbase[i] = a;
				tparamters[2 * i] = b - a;
				tparamters[2 * i + 1] = c - a;
				tnormal[i] = Point::cross(tparamters[2 * i], tparamters[2 * i + 1]).normalize();
			} break;
			case Element::CODE::SQUARE4:
//			case Element::CODE::SQUARE8:
			{
				const Point &a = surface->coordinates->datatarray()[enodes->at(0)];
				const Point &b = surface->coordinates->datatarray()[enodes->at(1)];
				const Point &c = surface->coordinates->datatarray()[enodes->at(2)];
				const Point &d = surface->coordinates->datatarray()[enodes->at(3)];
				Point center = (a + b + c + d) / 4;
				tnormal[i] = Point::cross(c - a, d - b).normalize();
				Point plane[4] = {
						a - tnormal[i] * (tnormal[i] * (a - center)),
						b - tnormal[i] * (tnormal[i] * (b - center)),
						c - tnormal[i] * (tnormal[i] * (c - center)),
						d - tnormal[i] * (tnormal[i] * (d - center))
				};
				Point pp[4] = { plane[1] - plane[0], plane[2] - plane[1], plane[3] - plane[2], plane[0] - plane[3] };
				for (int v = 0; v < 4; ++v) {
					double _s, _t;
					plane[(v + 2) % 4].getBarycentric(plane[v], pp[v], -pp[(v + 3) % 4], _s, _t);
					if (0 <= _s && _s <= 1 && 0 <= _t && _t <= 1) {
						tbase[i] = plane[v];
						tparamters[2 * i] = pp[v];
						tparamters[2 * i + 1] = -pp[(v + 3) % 4];
						break;
					}
				}
			} break;
			default:
				eslog::internalFailure("unknown or not implemented surface element.\n");
			}
		}

		normal[t].swap(tnormal);
		parameters[t].swap(tparamters);
		base[t].swap(tbase);
	}

	surface->normal = new serializededata<esint, Point>(1, normal);
	surface->parameters = new serializededata<esint, Point>(2, parameters);
	surface->base = new serializededata<esint, Point>(1, base);

	DebugOutput::warpedNormals("surface.planes", 1, 1);
	eslog::checkpointln("MESH: WARPED SURFACE NORMALS COMMPUTED");
	profiler::syncend("compute_warped_surface_normals");
}

void exchangeContactHalo()
{
	_Point<float> box[2] = {
			_Point<float>( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max(),  std::numeric_limits<float>::max()),
			_Point<float>(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max())
	};
	auto enodes = info::mesh->surface->enodes->begin();
	for (esint e = 0; e < info::mesh->surface->size; ++e, ++enodes) {
		if (info::mesh->surface->contact->datatarray()[e].gap > 0) {
			for (auto n = enodes->begin(); n != enodes->end(); ++n) {
				const Point &p = info::mesh->surface->coordinates->datatarray()[*n];
				box[0].x = std::min((float)p.x - info::mesh->surface->contact->datatarray()[e].gap, box[0].x);
				box[0].y = std::min((float)p.y - info::mesh->surface->contact->datatarray()[e].gap, box[0].y);
				box[0].z = std::min((float)p.z - info::mesh->surface->contact->datatarray()[e].gap, box[0].z);
				box[1].x = std::max((float)p.x + info::mesh->surface->contact->datatarray()[e].gap, box[1].x);
				box[1].y = std::max((float)p.y + info::mesh->surface->contact->datatarray()[e].gap, box[1].y);
				box[1].z = std::max((float)p.z + info::mesh->surface->contact->datatarray()[e].gap, box[1].z);
			}
		}
	}

	std::vector<_Point<float> > boxes(2 * info::mpi::size);
	Communication::allGather(box, boxes.data(), 6, MPI_FLOAT);

	auto areIntersected = [&] (_Point<float> *block1, _Point<float> *block2) {
		return	!
				(block1[1].x < block2[0].x || block2[1].x < block1[0].x) ||
				(block1[1].y < block2[0].y || block2[1].y < block1[0].y) ||
				(block1[1].z < block2[0].z || block2[1].z < block1[0].z);
	};

	info::mesh->contacts->neighbors.clear();
	for (int r = 0; r < info::mpi::size; r++) {
		if (r != info::mpi::rank) {
			if (areIntersected(box, boxes.data() + 2 * r)) {
				info::mesh->contacts->neighbors.push_back(r);
			}
		}
	}

	const auto &coordinates = info::mesh->surface->coordinates->datatarray();
	std::vector<std::vector<esint> > sBuffer(info::mesh->contacts->neighbors.size()), rBuffer(info::mesh->contacts->neighbors.size());
	for (size_t n = 0; n < info::mesh->contacts->neighbors.size(); ++n) {
		std::vector<esint> esend, nsend;
		const auto &min = boxes[2 * info::mesh->contacts->neighbors[n]];
		const auto &max = boxes[2 * info::mesh->contacts->neighbors[n] + 1];
		auto enodes = info::mesh->surface->enodes->begin();
		for (size_t e = 0; e < info::mesh->surface->enodes->structures(); ++e, ++enodes) {
			for (auto ec = enodes->begin(); ec != enodes->end(); ++ec) {
				if (
						min.x <= coordinates[*ec].x && coordinates[*ec].x <= max.x &&
						min.y <= coordinates[*ec].y && coordinates[*ec].y <= max.y &&
						min.z <= coordinates[*ec].z && coordinates[*ec].z <= max.z) {

					esend.push_back(e);
					for (auto c = enodes->begin(); c != enodes->end(); ++c) {
						nsend.push_back(*c);
					}
					break;
				}
			}
		}
		size_t ssize = 0, nsize = nsend.size();
		ssize += 1 + 3 * esend.size(); // fID, body, epointer size
		ssize += 1 + esend.size() + nsend.size(); // enodes size
		ssize += 1 + esend.size() * (4 * sizeof(Point) / sizeof(esint)); // plane
		utils::sortAndRemoveDuplicates(nsend);
		ssize += 1 + nsend.size() * (1 + sizeof(Point) / sizeof(esint)); // ids + coordinates size

		sBuffer[n].reserve(ssize);

		// send offset and bodies
		sBuffer[n].push_back(esend.size());
		for (size_t e = 0; e < esend.size(); ++e) {
			sBuffer[n].push_back(esend[e]);
		}
		for (size_t e = 0; e < esend.size(); ++e) {
			sBuffer[n].push_back(info::mesh->surface->body->datatarray()[esend[e]]);
		}
		for (size_t e = 0; e < esend.size(); ++e) {
			sBuffer[n].push_back(static_cast<int>(info::mesh->surface->epointers->datatarray()[esend[e]]->code));
		}

		// send enodes in target offsets
		sBuffer[n].push_back(esend.size() + nsize);
		enodes = info::mesh->surface->enodes->begin();
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
			const auto &nn = info::mesh->surface->normal->datatarray()[esend[e]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nn), reinterpret_cast<const esint*>(&nn) + sizeof(nn) / sizeof(esint));
		}
		for (size_t e = 0; e < esend.size(); ++e) {
			const auto &nn = info::mesh->surface->parameters->datatarray()[2 * esend[e]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nn), reinterpret_cast<const esint*>(&nn) + 2 * sizeof(nn) / sizeof(esint));
		}
		for (size_t e = 0; e < esend.size(); ++e) {
			const auto &nc = info::mesh->surface->base->datatarray()[esend[e]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nc), reinterpret_cast<const esint*>(&nc) + sizeof(nc) / sizeof(esint));
		}

		// send global ids
		sBuffer[n].push_back(nsend.size());
		for (size_t c = 0; c < nsend.size(); ++c) {
			sBuffer[n].push_back(info::mesh->surface->nIDs->datatarray()[nsend[c]]);
		}

		// send coordinates
		for (size_t c = 0; c < nsend.size(); ++c) {
			const auto &p = info::mesh->surface->coordinates->datatarray()[nsend[c]];
			sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&p), reinterpret_cast<const esint*>(&p) + sizeof(p) / sizeof(esint));
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->contacts->neighbors)) {
		eslog::internalFailure("cannot exchange contact halo.\n");
	}

	info::mesh->contacts->surfaces.resize(info::mesh->contacts->neighbors.size());
	for (size_t n = 0, i = 0; n < info::mesh->contacts->neighbors.size(); ++n, i = 0) {
		info::mesh->contacts->surfaces[n] = new SurfaceStore();
		// receive parents and bodies
		info::mesh->contacts->surfaces[n]->fID = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
		memcpy(info::mesh->contacts->surfaces[n]->fID->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(esint) * rBuffer[n][i]);
		info::mesh->contacts->surfaces[n]->body = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
		memcpy(info::mesh->contacts->surfaces[n]->body->datatarray().data(), rBuffer[n].data() + i + 1 + rBuffer[n][i], sizeof(esint) * rBuffer[n][i]);
		info::mesh->contacts->surfaces[n]->epointers = new serializededata<esint, Element*>(1, tarray<Element*>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
		for (esint e = 0; e < rBuffer[n][i]; ++e) {
			info::mesh->contacts->surfaces[n]->epointers->datatarray()[e] = &Mesh::edata[rBuffer[n][2 * rBuffer[n][i] + e + i + 1]];
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
		serializededata<esint, esint>::balance(edist, enodes, &info::mesh->contacts->surfaces[n]->fID->datatarray().distribution());
		info::mesh->contacts->surfaces[n]->enodes = new serializededata<esint, esint>(edist, enodes);

		// receive warped normals
		size = rBuffer[n][i];
		info::mesh->contacts->surfaces[n]->normal = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->surfaces[n]->normal->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i + 1), sizeof(Point) * size);
		i += 1 + size * (sizeof(Point) / sizeof(esint));

		info::mesh->contacts->surfaces[n]->parameters = new serializededata<esint, Point>(2, tarray<Point>(info::env::OMP_NUM_THREADS, 2 * size));
		memcpy(info::mesh->contacts->surfaces[n]->parameters->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i), 2 * sizeof(Point) * size);
		i += size * (2 * sizeof(Point) / sizeof(esint));

		info::mesh->contacts->surfaces[n]->base = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->surfaces[n]->base->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i), sizeof(Point) * size);
		i += size * (sizeof(Point) / sizeof(esint));

		// receive global ids
		size = rBuffer[n][i];
		info::mesh->contacts->surfaces[n]->nIDs = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->surfaces[n]->nIDs->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(esint) * rBuffer[n][i]);
		i += 1 + rBuffer[n][i];

		// receive coordinates
		info::mesh->contacts->surfaces[n]->coordinates = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
		memcpy(info::mesh->contacts->surfaces[n]->coordinates->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i), sizeof(Point) * size);
	}
	info::mesh->contacts->neighborsWithMe = info::mesh->contacts->neighbors;
	info::mesh->contacts->neighborsWithMe.push_back(info::mpi::rank);
	info::mesh->contacts->surfaces.push_back(info::mesh->surface);

	eslog::checkpointln("MESH: CLOSE BOUNDARY EXCHANGED");
}

void findCloseElements()
{
	// checking
	// 1. distance to plane defined by normal and center
	// 2. test if any point is in coarse element defined by: base + s * parameter.u + t * parameter.v

	ContactStore *contact = info::mesh->contacts;

	// neighbors, local
	std::vector<Point> nstart, nend, lstart, lend;
	std::vector<esint> neigh, offset = { 0 };

	for (size_t n = 0; n < contact->neighborsWithMe.size() - 1; ++n) {
		offset.push_back(offset.back() + contact->surfaces[n]->enodes->structures());
	}
	// push neighbors
	nstart.reserve(offset.back());
	nend.reserve(offset.back());
	neigh.reserve(offset.back());
	for (size_t r = 0, offset = 0; r < contact->neighborsWithMe.size() - 1; ++r) {
		for (auto e = contact->surfaces[r]->enodes->begin(); e != contact->surfaces[r]->enodes->end(); ++e, ++offset) {
			nstart.push_back(contact->surfaces[r]->coordinates->datatarray()[e->front()]);
			nend.push_back(nstart.back());
			for (auto n = e->begin() + 1; n != e->end(); ++n) {
				contact->surfaces[r]->coordinates->datatarray()[*n].minmax(nstart.back(), nend.back());
			}
		}
		neigh.insert(neigh.end(), contact->surfaces[r]->enodes->structures(), r);
	}
	// push local
	lstart.reserve(contact->surfaces.back()->enodes->structures());
	lend.reserve(contact->surfaces.back()->enodes->structures());
	for (auto e = contact->surfaces.back()->enodes->begin(); e != contact->surfaces.back()->enodes->end(); ++e) {
		lstart.push_back(contact->surfaces.back()->coordinates->datatarray()[e->front()]);
		lend.push_back(lstart.back());
		for (auto n = e->begin() + 1; n != e->end(); ++n) {
			contact->surfaces.back()->coordinates->datatarray()[*n].minmax(lstart.back(), lend.back());
		}
	}

	IntervalTree ntree(nstart, nend);
	IntervalTree ltree(lstart, lend);

	std::unordered_map<esint, std::unordered_set<esint> > bodyPairs;
	std::vector<esint> dist = { 0 }, data;
	dist.reserve(contact->surfaces.back()->size + 1);
	data.reserve(offset.back() + contact->surfaces.back()->size);
	std::vector<esint> nintervals, lintervals;

	for (esint face = 0; face < contact->surfaces.back()->size; ++face) {
		bool self_contact = contact->surfaces.back()->contact->datatarray()[face].self_contact;
		double gap = contact->surfaces.back()->contact->datatarray()[face].gap;
		double angle = -std::cos(M_PI * contact->surfaces.back()->contact->datatarray()[face].angle / 180);

		int body = contact->surfaces.back()->body->datatarray()[face];
		Point normal = contact->surfaces.back()->normal->datatarray()[face];
		Point base = contact->surfaces.back()->base->datatarray()[face];
		Point u = contact->surfaces.back()->parameters->datatarray()[2 * face];
		Point v = contact->surfaces.back()->parameters->datatarray()[2 * face + 1];
		bool istriangle =
				(contact->surfaces.back()->epointers->datatarray()[face]->code == Element::CODE::TRIANGLE3) ||
				(contact->surfaces.back()->epointers->datatarray()[face]->code == Element::CODE::TRIANGLE6);

		auto cohen_sutherland = [&] (int n, esint index) {
			auto getCode = [&] (const double &d, const double &s, const double &t) {
				int code = 0;
				code |= (s < 0) ? 1 : 0;
				code |= (t < 0) ? 2 : 0;
				if (istriangle) {
					code |= (1 < t + s) ? 4 : 0;
				} else {
					code |= (1 < s) ? 4 : 0;
					code |= (1 < t) ? 8 : 0;
				}
				code |= d < -gap ? 16 : 0;
				code |= gap < d  ? 32 : 0;
				return code;
			};

			int code = 63;
			auto nodes = contact->surfaces[n]->enodes->begin() + index;
			for (auto nn = nodes->begin(); nn != nodes->end(); ++nn) {
				Point p = contact->surfaces[n]->coordinates->datatarray()[*nn];
				double s, t, d = normal * (p - base);
				p -= normal * d;
				p.getBarycentric(base, u, v, s, t);
				code &= getCode(d, s, t);
			}

			return code == 0;
		};

		Point gapdirection = normal;
		gapdirection.abs();
		nintervals.clear();
		lintervals.clear();
		ntree.traverse(1, lstart[face] - gapdirection * gap, lend[face] + gapdirection * gap, nintervals);
		ltree.traverse(1, lstart[face] - gapdirection * gap, lend[face] + gapdirection * gap, lintervals);

		for (size_t i = 0; i < nintervals.size(); ++i) {
			for (auto pp = ntree.permutation.cbegin() + ntree.begin(nintervals[i]); pp != ntree.permutation.cbegin() + ntree.end(nintervals[i]); ++pp) {
				int n = neigh[*pp];
				esint index = *pp - offset[n];
				int nbody = contact->surfaces[n]->body->datatarray()[index];
				const Point &nnormal = contact->surfaces[n]->normal->datatarray()[index];

				if ((self_contact || nbody != body) && nnormal * normal <= angle && cohen_sutherland(n, index)) {
					bodyPairs[std::min(body, nbody)].insert(std::max(body, nbody));
					data.push_back(n);
					data.push_back(index);
				}
			}
		}

		for (size_t i = 0; i < lintervals.size(); ++i) {
			for (auto pp = ltree.permutation.cbegin() + ltree.begin(lintervals[i]); pp != ltree.permutation.cbegin() + ltree.end(lintervals[i]); ++pp) {
				int nbody = contact->surfaces.back()->body->datatarray()[*pp];
				const Point &nnormal = contact->surfaces.back()->normal->datatarray()[*pp];
//				printf("%d->%d=%d-%d%d%d\n", face, *pp, *pp < face, (self_contact || nbody != body), nnormal * normal <= angle, cohen_sutherland(contact->surfaces.size() - 1, *pp));
				if ((self_contact || nbody != body) && nnormal * normal <= angle && cohen_sutherland(contact->surfaces.size() - 1, *pp)) {
					bodyPairs[std::min(body, nbody)].insert(std::max(body, nbody));
					data.push_back(contact->surfaces.size() - 1);
					data.push_back(*pp);
				}
			}
		}

		dist.push_back(data.size());
	}

	contact->pairs = new serializededata<esint, esint>(dist, data);

	std::vector<std::pair<esint, esint> > interfaces;
	for (auto i = bodyPairs.begin(); i != bodyPairs.end(); ++i) {
		for (auto j = i->second.begin(); j != i->second.end(); ++j) {
			interfaces.push_back(std::make_pair(i->first, *j));
		}
	}
	std::sort(interfaces.begin(), interfaces.end());
	Communication::uniqueAllGatherUnknownSize(interfaces);
	for (auto i = interfaces.begin(); i != interfaces.end(); ++i) {
		contact->interfaces.push_back(Interface(i->first, i->second));
	}

	eslog::checkpointln("MESH: CLOSE ELEMENTS FOUND");
}

void triangulate(std::vector<Point> &face, std::vector<Triangle> &triangles)
{
	triangles.clear();
	auto prev= [&] (const int &i) -> int {
		return (i + face.size() - 1) % face.size();
	};
	auto next = [&] (const int &i) -> int {
		return (i + 1) % face.size();
	};

	if (face.size() == 3) {
		triangles.push_back(Triangle{ face, 0, 1, 2 });
		return;
	}
	if (face.size() == 4) {
		for (int v = 0, size = face.size(); v < size; ++v) {
			Point v1 = face[v] - face[prev(v)], v2 = face[next(v)] - face[prev(v)];
			if (Point::cross2d(v1.normalize(), v2.normalize()) < 0) {
				triangles.push_back(Triangle{ face, v, next(v), next(next(v)) });
				triangles.push_back(Triangle{ face, next(next(v)), prev(v), v });
				return;
			}
		}
		if ((face[0] - face[2]).length() < (face[1] - face[3]).length()) {
			triangles.push_back(Triangle{ face, 0, 1, 2 });
			triangles.push_back(Triangle{ face, 0, 2, 3 });
		} else {
			triangles.push_back(Triangle{ face, 0, 1, 3 });
			triangles.push_back(Triangle{ face, 1, 2, 3 });
		}
		return;
	}
	eslog::internalFailure("cannot compute triangles from a polygon.\n");
}

void clip(const Point &base, const std::vector<Triangle> &triangles, const std::vector<Point> &polygon, const double &gap, std::vector<Triangle> &output)
{
	double eps = 1e-10;
	enum status { in, out, on, processed };
	std::vector<Triangle> res;

	auto crossPoint = [] (const Point &p0, const Point &p1, const Point &q0, const Point &q1) {
		Point u = p1 - p0, v = q1 - q0, w = p0 - q0;
		double t = Point::cross2d(u, w) / Point::cross2d(u, v);
		return q0 + v * t;
	};

	for (size_t i = 0; i < triangles.size(); ++i) {
		std::vector<Point> in(polygon.rbegin(), polygon.rend()), nextin;
		auto curr = [&] (size_t i) { return (i                ) % in.size(); };
		auto next = [&] (size_t i) { return (i + 1            ) % in.size(); };
		auto prev = [&] (size_t i) { return (i - 1 + in.size()) % in.size(); };

		Point p[3] = { triangles[i].p[0], triangles[i].p[1], triangles[i].p[2] };
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
					nextin.push_back(crossPoint(p[c], p[(c + 1) % 3], in[prev(j)], in[j]));
					status[prev(j)] = status::processed;
					for (size_t k = 0; k < in.size(); ++k) {
						if (status[curr(j + k)] == status::in || status[curr(j + k)] == status::on) {
							nextin.push_back(in[curr(j + k)]);
							status[curr(j + k)] = status::processed;
						} else {
							nextin.push_back(crossPoint(p[c], p[(c + 1) % 3], in[prev(j + k)], in[curr(j + k)]));
							status[curr(j + k)] = status::processed;
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
			std::vector<Point> _res;
			for (size_t j = 0; j < in.size(); ++j) {
				if (std::fabs(in[prev(j)].x - in[j].x) > 1e-6 || std::fabs(in[prev(j)].y - in[j].y) > 1e-6) {
					_res.push_back(in[j]);
				}
			}
			if (3 <= _res.size()) {
				for (esint j = 2, size = _res.size(); j < size; j++) { // Delaunay triangulation?
					res.push_back(Triangle{ _res, 0, j - 1, j });
				}
			}
		}
	}

	// remove triangles that are above the limit
	auto max = [&] (const double &z, const Point &p0, const Point &p1) { return p0 + (p1 - p0) * ((z - p0.z) / (p1.z - p0.z)); };

	output.clear();
	for (size_t t = 0; t < res.size(); ++t) {
		bool below[3] = { res[t].p[0].z + gap < base.z, res[t].p[1].z + gap < base.z, res[t].p[2].z + gap < base.z };

		if (!below[0] && !below[1] && !below[2]) {
			output.push_back(res[t]);
		} else if (!below[0] || !below[1] || !below[2]) {
			int i = 0;
			while (below[0] || !below[2]) { std::rotate(below, below + 1, below + 3); ++i; }
			if (below[1]) {
				output.push_back(Triangle(res[t].p[i], max(base.z - gap, res[t].p[i], res[t].p[(i + 1) % 3]), max(base.z - gap, res[t].p[i], res[t].p[(i + 2) % 3])));
			} else {
				Point c0 = max(base.z - gap, res[t].p[i], res[t].p[(i + 2) % 3]);
				Point c1 = max(base.z - gap, res[t].p[(i + 1) % 3], res[t].p[(i + 2) % 3]);
				output.push_back(Triangle(res[t].p[i], res[t].p[(i + 1) % 3], c0));
				output.push_back(Triangle(res[t].p[(i + 1) % 3], c1, c0));
			}
		}
	}
	output.swap(res);
	output.clear();
	for (size_t t = 0; t < res.size(); ++t) {
		bool above[3] = { base.z < res[t].p[0].z - gap, base.z < res[t].p[1].z - gap, base.z < res[t].p[2].z - gap };

		if (!above[0] && !above[1] && !above[2]) {
			output.push_back(res[t]);
		} else if (!above[0] || !above[1] || !above[2]) {
			int i = 0;
			while (above[0] || !above[2]) { std::rotate(above, above + 1, above + 3); ++i; }
			if (above[1]) {
				output.push_back(Triangle(res[t].p[i], max(base.z + gap, res[t].p[i], res[t].p[(i + 1) % 3]), max(base.z + gap, res[t].p[i], res[t].p[(i + 2) % 3])));
			} else {
				Point c0 = max(base.z + gap, res[t].p[i], res[t].p[(i + 2) % 3]);
				Point c1 = max(base.z + gap, res[t].p[(i + 1) % 3], res[t].p[(i + 2) % 3]);
				output.push_back(Triangle(res[t].p[i], res[t].p[(i + 1) % 3], c0));
				output.push_back(Triangle(res[t].p[(i + 1) % 3], c1, c0));
			}
		}
	}
}

void computeContactInterface()
{
	profiler::syncstart("compute_contact_interface");

	const std::vector<SurfaceStore*> &surfaces = info::mesh->contacts->surfaces;

	Point axis;
	double cos, sin;
	auto setRotation = [&] (esint e) {
		const Point &normal = info::mesh->surface->normal->datatarray()[e];
		if (normal.z < -0.999) {
			axis.x = 1;
			sin = axis.y = axis.z = 0;
			cos = -1;
		} else if (normal.z < 0.999) {
			axis = Point::cross(normal, Point(0, 0, 1)).normalize();
			double angle = std::acos(normal.z);
			sin = std::sin(angle);
			cos = std::cos(angle);
		} else {
			sin = axis.x = axis.y = 0;
			cos = axis.z = 1;
		}
	};

	std::vector<Point> plane(20), projected(20);
	std::vector<Triangle> planeTriangles, intersection;
	auto setPolygon = [&] (std::vector<Point> &polygon, esint neigh, esint offset) {
		polygon.clear();
		auto nodes = surfaces[neigh]->enodes->cbegin() + offset;
		const auto &epointer = surfaces[neigh]->epointers->datatarray();
		for (auto n = epointer[offset]->polygon->begin(); n != epointer[offset]->polygon->end(); ++n) {
			polygon.push_back(surfaces[neigh]->coordinates->datatarray()[nodes->at(*n)]);
			polygon.back().rodrigues(axis, cos, sin);
		}
	};

	std::vector<Triangle> triangles;

	std::vector<esint> sdist = { 0 }, ddist = { 0 }, pdist = { 0 };
	std::vector<SparseSegment> sparse;
	std::vector<DenseSegment> dense;
	std::vector<Point2D> planeCoordinates;

	auto pairs = info::mesh->contacts->pairs->cbegin();
	for (esint e = 0; e < info::mesh->surface->size; ++e, ++pairs) {
		if (pairs->size() == 0) {
			continue;
		}

		setRotation(e);

		double earea = 0;
		setPolygon(plane, info::mesh->contacts->neighbors.size(), e);
		Point base = info::mesh->surface->base->datatarray()[e];
		base.rodrigues(axis, cos, sin);
		triangulate(plane, planeTriangles);
		for (size_t t = 0; t < planeTriangles.size(); ++t) {
			earea += planeTriangles[t].area();
		}

		std::unordered_map<esint, double> insertedBodies;
		for (auto other = pairs->begin(); other != pairs->end(); ++other) {
			esint neigh = *other++;
			esint offset = *other;

			setPolygon(projected, neigh, *other);
			clip(base, planeTriangles, projected, info::mesh->surface->contact->datatarray()[e].gap, intersection);
			if (intersection.size()) {
				if (insertedBodies.empty()) {
					sparse.push_back(SparseSegment(surfaces.back()->body->datatarray()[e], e, planeCoordinates.size(), triangles.size(), dense.size()));
					planeCoordinates.insert(planeCoordinates.end(), plane.begin(), plane.end());
				}
				if (insertedBodies.count(surfaces[neigh]->body->datatarray()[offset]) == 0) {
					insertedBodies[surfaces[neigh]->body->datatarray()[offset]] = 0;
				}

				++sparse.back().denseSegmentEnd;
				dense.push_back(DenseSegment(neigh, surfaces[neigh]->body->datatarray()[offset], offset, planeCoordinates.size(), planeCoordinates.size() + projected.size()));
				planeCoordinates.insert(planeCoordinates.end(), projected.begin(), projected.end());
				for (size_t i = 0; i < intersection.size(); i++) {
					planeCoordinates.insert(planeCoordinates.end(), intersection[i].p, intersection[i].p + 3);
				}
				dense.back().triangles += intersection.size();
			}
			for (size_t i = 0; i < intersection.size(); i++) {
				double area = intersection[i].area();
				insertedBodies[surfaces[neigh]->body->datatarray()[offset]] += area / earea;
				intersection[i].rotate(axis, cos, -sin);
				triangles.push_back(intersection[i]);
			}
		}
		for (auto ib = insertedBodies.begin(); ib != insertedBodies.end(); ++ib) {
			if (ib->second < MIN_SLAVE_COVER_RATIO) {
				for (esint d = sparse.back().denseSegmentBegin; d < sparse.back().denseSegmentEnd; ++d) {
					if (surfaces[dense[d].neigh]->body->datatarray()[dense[d].element] == ib->first) {
						dense[d].skip = true;
					}
				}
			}
		}
		if (insertedBodies.size()) {
			sdist.push_back(sparse.size());
			ddist.push_back(dense.size());
			pdist.push_back(planeCoordinates.size());
		}
	}

	info::mesh->contacts->sparseSide = new serializededata<esint, SparseSegment>(sdist, sparse);
	info::mesh->contacts->denseSide = new serializededata<esint, DenseSegment>(ddist, dense);
	info::mesh->contacts->planeCoordinates = new serializededata<esint, Point2D>(pdist, planeCoordinates);
	info::mesh->contacts->intersections = new serializededata<esint, Triangle>(1, triangles);

	profiler::syncend("compute_contact_interface");
	eslog::checkpointln("MESH: CONTACT INTERFACE COMPUTED");
}

void arrangeContactInterfaces()
{
	profiler::syncstart("arrange_contact_interface");

	const std::vector<SurfaceStore*> &surfaces = info::mesh->contacts->surfaces;

	struct istats { double area; esint faces, triangles; bool skip; istats(): area(0), faces(0), triangles(0), skip(true) {} };
	std::unordered_map<esint, std::unordered_map<esint, istats> > istats;

	auto *sside = info::mesh->contacts->sparseSide;
	auto *dside = info::mesh->contacts->denseSide;
	auto &coors = info::mesh->contacts->planeCoordinates->datatarray();
	for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
		bool add = true;
		for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
			if (d->skip) {
				continue;
			}

			if (add) {
				add = false;
				istats[s->body][d->body].faces += 1;
			}
			istats[s->body][d->body].triangles += d->triangles;

			for (esint t = 0; t < d->triangles; ++t) {
				istats[s->body][d->body].area += Triangle::area(coors.data() + d->triangleOffset + 3 * t);
			}
		}
	}

	std::vector<Interface> interfaces = info::mesh->contacts->interfaces;
	auto comp = [] (const Interface &i, const Interface &j) {
		if (i.from.body == j.from.body) {
			return i.to.body < j.to.body;
		}
		return i.from.body < j.from.body;
	};
	auto size = interfaces.size();
	std::vector<double> area(2 * size);
	std::vector<esint> faces(2 * size);
	std::vector<esint> triasum(2 * size), triaoffset(2 * size);
	for (auto i = istats.begin(); i != istats.end(); ++i) {
		for (auto j = i->second.begin(); j != i->second.end(); ++j) {
			auto it = std::lower_bound(interfaces.begin(), interfaces.end(), Interface(i->first, j->first), comp);
			if (it == interfaces.end() || it->from.body != i->first || it->to.body != j->first) {
				it = std::lower_bound(interfaces.begin(), interfaces.end(), Interface(j->first, i->first), comp);
				auto offset = it - interfaces.begin();
				area[offset + size] += j->second.area;
				faces[offset + size] += j->second.faces;
				triasum[offset + size] += j->second.triangles;
			} else {
				auto offset = it - interfaces.begin();
				area[offset] += j->second.area;
				faces[offset] += j->second.faces;
				triasum[offset] += j->second.triangles;
			}
		}
	}

	Communication::allReduce(area, Communication::OP::SUM);
	Communication::allReduce(faces, Communication::OP::SUM);
	std::vector<esint> triasize = triaoffset = triasum;
	Communication::exscan(triasum, triaoffset);
	info::mesh->contacts->interfaces.clear();
	for (size_t i = 0; i < interfaces.size(); ++i) {
		if (faces[i] || faces[i + size]) {
			info::mesh->contacts->interfaces.push_back(Interface(interfaces[i].from.body, interfaces[i].to.body));
			info::mesh->contacts->interfaces.back().from.area = area[i];
			info::mesh->contacts->interfaces.back().from.faces = faces[i];
			info::mesh->contacts->interfaces.back().from.triangleOffset = triaoffset[i];
			info::mesh->contacts->interfaces.back().from.triangleSize = triasize[i];
			info::mesh->contacts->interfaces.back().from.triangleTotalSize = triasum[i];
			info::mesh->contacts->interfaces.back().to.area = area[i + size];
			info::mesh->contacts->interfaces.back().to.faces = faces[i + size];
			info::mesh->contacts->interfaces.back().to.triangleOffset = triaoffset[i + size];
			info::mesh->contacts->interfaces.back().to.triangleSize = triasize[i + size];
			info::mesh->contacts->interfaces.back().to.triangleTotalSize = triasum[i + size];
			info::mesh->contacts->interfaces.back().setOrientation();
		}
	}
	std::sort(info::mesh->contacts->interfaces.begin(), info::mesh->contacts->interfaces.end(), comp);

	std::vector<std::string> bnames(info::mesh->bodies->totalSize);
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r) {
		for (size_t b = 0; b < info::mesh->elementsRegions[r]->bodies.size(); ++b) {
			if (bnames[info::mesh->elementsRegions[r]->bodies[b]].size()) {
				bnames[info::mesh->elementsRegions[r]->bodies[b]] += "_";
			}
			bnames[info::mesh->elementsRegions[r]->bodies[b]] += info::mesh->elementsRegions[r]->name;
		}
	}

	esint myrank = surfaces.size() - 1;
	std::vector<int> assigned(info::mesh->contacts->interfaces.size()), current;
	std::vector<std::vector<esint> > sBuffer(info::mesh->contacts->neighbors.size(), std::vector<esint>(info::mesh->contacts->interfaces.size() + 1)), rBuffer(info::mesh->contacts->neighbors.size());

	auto preprocess = [&] (size_t index, esint from, esint to) {
		istats[from][to].skip = false;
		for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
			if (s->body != from) {
				continue;
			}
			for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
				if (d->body == to && d->neigh < myrank) {
					++sBuffer[d->neigh][index];
					sBuffer[d->neigh].push_back(surfaces[d->neigh]->fID->datatarray()[d->element]);
				}
			}
		}
	};

	auto create = [&] (const std::string &name, size_t index, esint from, esint to) {
		info::mesh->contactInterfaces.push_back(new ContactInterfaceStore("CONTACT-" + name + "-" + bnames[from] + "-" + bnames[to], index));

		std::vector<esint> dist = { 0 }, data;
		std::vector<Element*> epointer;

		auto push = [&] (esint e) {
			auto face = surfaces.back()->enodes->begin() + e;
			epointer.push_back(surfaces.back()->epointers->datatarray()[e]);
			for (auto n = face->begin(); n != face->end(); ++n) {
				data.push_back(surfaces.back()->nodes->datatarray()[*n]);
			}
			dist.push_back(data.size());
		};

		std::vector<esint> dense;
		for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
			if (s->body != from) {
				continue;
			}
			bool addsparse = true;
			for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
				if (d->body == to) {
					if (addsparse) {
						addsparse = false;
						push(s->element);
					}
					if (d->neigh == myrank) {
						dense.push_back(d->element);
					}
				}
			}
		}
		for (size_t n = 0; n < rBuffer.size(); ++n) {
			for (esint i = rBuffer[n][index]; i < rBuffer[n][index + 1]; ++i) {
				dense.push_back(rBuffer[n][i]);
			}
		}
		utils::sortAndRemoveDuplicates(dense);
		for (size_t i = 0; i < dense.size(); ++i) {
			push(dense[i]);
		}

		info::mesh->contactInterfaces.back()->originalDimension = 2;
		info::mesh->contactInterfaces.back()->dimension = 2;
		info::mesh->contactInterfaces.back()->procNodes = new serializededata<esint, esint>(dist, data);
		info::mesh->contactInterfaces.back()->epointers = new serializededata<esint, Element*>(1, epointer);
	};

	for (int insert = 0; insert <= 1; ++insert) {
		for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
			switch (it->second.detection) {
			case ContactInterfaceConfiguration::DETECTION::ALL_BODIES:
				for (size_t c = 0; c < info::mesh->contacts->interfaces.size(); ++c) {
					if (assigned[c] == insert) {
						if (insert) {
							it->second.found_interfaces.push_back(info::mesh->contactInterfaces.size());
							create(it->first, c, info::mesh->contacts->interfaces[c].from.body, info::mesh->contacts->interfaces[c].to.body);
						} else {
							preprocess(c, info::mesh->contacts->interfaces[c].from.body, info::mesh->contacts->interfaces[c].to.body);
						}
						assigned[c] = insert + 1;
					}
				}
				break;
			case ContactInterfaceConfiguration::DETECTION::BODY_LIST:
				current.clear();
				for (size_t b = 0; b < it->second.body_list.size(); ++b) {
					const ElementsRegionStore *region = info::mesh->eregion(it->second.body_list[b]);
					current.insert(current.end(), region->bodies.begin(), region->bodies.end());
				}
				utils::sortAndRemoveDuplicates(current);
				for (size_t c = 0; c < info::mesh->contacts->interfaces.size(); ++c) {
					if (
							assigned[c] == insert &&
							std::binary_search(current.begin(), current.end(), info::mesh->contacts->interfaces[c].from.body) &&
							std::binary_search(current.begin(), current.end(), info::mesh->contacts->interfaces[c].to.body)) {

						if (insert) {
							create(it->first, c, info::mesh->contacts->interfaces[c].from.body, info::mesh->contacts->interfaces[c].to.body);
						} else {
							preprocess(c, info::mesh->contacts->interfaces[c].from.body, info::mesh->contacts->interfaces[c].to.body);
						}
						assigned[c] = insert + 1;
					}
				}
				break;
			case ContactInterfaceConfiguration::DETECTION::CONTACT_PAIR:
				break;
			}
		}

		if (insert == 0) {
			for (size_t n = 0; n < sBuffer.size(); ++n) {
				esint sum = info::mesh->contacts->interfaces.size() + 1;
				for (size_t i = 0; i <= info::mesh->contacts->interfaces.size(); ++i) {
					esint tmp = sBuffer[n][i];
					sBuffer[n][i] = sum;
					sum += tmp;
				}
			}
			if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->contacts->neighbors)) {
				eslog::internalFailure("cannot exchange contact faces to neighbors.\n");
			}
		}
	}

	for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
		for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
			d->skip = d->skip | istats[surfaces.back()->body->datatarray()[s->element]][surfaces[d->neigh]->body->datatarray()[d->element]].skip;
		}
	}

	DebugOutput::contact(1, 1);
	profiler::syncend("arrange_contact_interface");
	eslog::checkpointln("MESH: CONTACT INTERFACE ARRANGED");
}

}
}
