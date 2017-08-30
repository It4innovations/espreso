
#include <fstream>

#include "mesh.h"
#include "mkl.h"

#include <numeric>

#include "../../basis/logging/logging.h"
#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../settings/evaluator.h"
#include "coordinates.h"
#include "material.h"
#include "region.h"

#include "../elements/point/node.h"
#include "../elements/point/dof.h"
#include "../elements/point/unknownpoint.h"
#include "../elements/line/line2.h"
#include "../elements/line/line3.h"
#include "../elements/line/unknownline.h"
#include "../elements/plane/square4.h"
#include "../elements/plane/square8.h"
#include "../elements/plane/triangle3.h"
#include "../elements/plane/triangle6.h"
#include "../elements/plane/unknownplane.h"

#include "../elements/element.h"
#include "elementtypes.h"

#include "metis.h"
#include "../../configuration/environment.h"
#include "../../configuration/material/coordinatesystem.h"

namespace espreso {

Mesh::Mesh(): _continuous(true), _elements(0)
{
	_coordinates = new Coordinates();
	_steps = 1;
	_bodies = { 0, 0 };
	_partPtrs = { 0, 0 };

	_regions.push_back(new Region(ElementType::ELEMENTS, _elements));
	_regions.back()->name = "ALL_ELEMENTS";
	_regions.push_back(new Region(ElementType::NODES, _nodes));
	_regions.back()->name = "ALL_NODES";
}

APIMesh::APIMesh(eslocal *l2g, size_t size)
: _l2g(l2g, l2g + size)
{
	_g2l = new std::vector<G2L>();
};

void Mesh::computeFixPoints(size_t number)
{
	if (_fixPoints.size() && _fixPoints[0].size() == number) {
		return;
	}

	_fixPoints.resize(parts());

	#pragma omp parallel for
	for  (size_t part = 0; part < parts(); part++) {
		size_t max = _partPtrs[part + 1] - _partPtrs[part];
		size_t points = std::min(number, max);
		std::vector<eslocal> fixPoints(points);
		std::vector<eslocal> eSubPartition = getPartition(_partPtrs[part], _partPtrs[part + 1], points);

		for (size_t j = 0; j < points; j++) {
			fixPoints[j] = getCentralNode(_partPtrs[part], _partPtrs[part + 1], eSubPartition, part, j);
		}
		std::sort(fixPoints.begin(), fixPoints.end());

		// Remove the same points
		Esutils::removeDuplicity(fixPoints);

		_fixPoints[part].reserve(fixPoints.size());
		for (size_t i = 0; i < fixPoints.size(); i++) {
			_fixPoints[part].push_back(_nodes[fixPoints[i]]);
		}
	}
}

static void checkMETISResult(eslocal result)
{
	switch (result) {
	case METIS_ERROR_INPUT:
		ESINFO(ERROR) << "An input for METIS procedure is incorrect.\n";
	case METIS_ERROR_MEMORY:
		ESINFO(ERROR) << "There is not enough memory for compute a partition.\n";
	case METIS_ERROR:
		ESINFO(ERROR) << "METIS fail computation.\n";
	}
}

static std::vector<eslocal> computePartition(std::vector<eslocal> &elements, std::vector<eslocal> &nodes, eslocal eSize, eslocal nSize, eslocal nCommon, eslocal parts)
{
	if (parts == 1) {
		return std::vector<eslocal> (elements.size(), 0);
	}
	// INPUTS
	eslocal options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_CONTIG] = 1;

	// OUTPUTS
	eslocal objval;
	std::vector<eslocal> ePartition(eSize), nPartition(nSize);

	checkMETISResult(METIS_PartMeshDual(
			&eSize,
			&nSize,
			elements.data(),
			nodes.data(),
			NULL,		// weights of nodes
			NULL,		// size of nodes
			&nCommon,
			&parts,
			NULL,		// weights of parts
			options,
			&objval,
			ePartition.data(),
			nPartition.data()));

	return ePartition;
}

static eslocal localIndex(const std::vector<eslocal> &nodeProjection, eslocal index)
{
	return std::lower_bound(nodeProjection.begin(), nodeProjection.end(), index) - nodeProjection.begin();
}

static void METISMeshRepresentation(
		const std::vector<Element*> &elements, size_t begin, size_t end,
		std::vector<eslocal> &nodeProjection, std::vector<eslocal> &metisElements, std::vector<eslocal> &metisNodes, eslocal &nCommon)
{
	for (size_t e = begin; e < end; e++) {
		for (size_t n = 0; n < elements[e]->coarseNodes(); n++) {
			nodeProjection.push_back(elements[e]->node(n));
		}
	}

	std::sort(nodeProjection.begin(), nodeProjection.end());
	Esutils::removeDuplicity(nodeProjection);

	// create array storing pointers to elements' nodes
	metisElements.reserve(end - begin + 1);
	metisElements.push_back(0);
	for (size_t i = begin; i < end; i++) {
		metisElements.push_back(metisElements.back() + elements[i]->coarseNodes());
	}

	// create array of nodes
	metisNodes.reserve(metisElements.back());
	// number of common nodes to be neighbor
	nCommon = 4;
	for (size_t i = begin; i < end; i++) {
		for (size_t j = 0; j < elements[i]->coarseNodes(); j++) {
			metisNodes.push_back(localIndex(nodeProjection, elements[i]->node(j)));
		}
		if (nCommon > elements[i]->nCommon()) {
			nCommon = elements[i]->nCommon();
		}
	}
}

static std::vector<eslocal> continuousReorder(std::vector<Element*> &elements, size_t begin, size_t end)
{
	eslocal nCommon;
	std::vector<eslocal> nodeProjection, metisElements, metisNodes;

	METISMeshRepresentation(elements, begin, end, nodeProjection, metisElements, metisNodes, nCommon);

	eslocal ne = end - begin, nn = nodeProjection.size(), numflag = 0, *xAdj, *adjncy;
	checkMETISResult(METIS_MeshToDual(&ne, &nn, metisElements.data(), metisNodes.data(), &nCommon, &numflag, &xAdj, &adjncy));

	std::vector<eslocal> partPtrs;
	std::vector<Element*> eCopy(elements.begin() + begin, elements.begin() + end);

	partPtrs.push_back(begin);
	size_t back = begin;
	for (size_t e = 0; e < eCopy.size(); e++) {
		if (eCopy[e] != NULL) {
			std::vector<eslocal> stack;
			elements[back++] = eCopy[e];
			eCopy[e] = NULL;
			stack.push_back(e);
			while (stack.size()) {
				eslocal current = stack.back();
				stack.pop_back();
				for (eslocal i = xAdj[current]; i < xAdj[current + 1]; i++) {
					if (eCopy[adjncy[i]] != NULL) {
						elements[back++] = eCopy[adjncy[i]];
						eCopy[adjncy[i]] = NULL;
						stack.push_back(adjncy[i]);
					}
				}
			}
			partPtrs.push_back(back);
		}
	}

	METIS_Free(xAdj);
	METIS_Free(adjncy);

	return partPtrs;
}

void Mesh::partitiateNoncontinuously(size_t parts, size_t noncontinuousParts)
{
	ESTEST(MANDATORY) << "Number of domains cannot be " << parts << (parts == 0 ? TEST_FAILED : TEST_PASSED);
	_continuousPartId.clear();
	_continuousPartId.resize(parts, 0);
	if (parts == 1 && this->parts() == 1) {
		_partPtrs = { 0, (eslocal)_elements.size() };
		mapElementsToDomains();
	} else {
		size_t blockSize = _elements.size() / noncontinuousParts;
		std::vector<eslocal> blocks;
		for (size_t b = 0; b < noncontinuousParts; b++) {
			std::vector<eslocal> reorder = continuousReorder(_elements, b * blockSize, b + 1 == noncontinuousParts ? _elements.size() : (b + 1) * blockSize);
			blocks.insert(blocks.end(), reorder.begin(), reorder.end() - 1);
		}
		blocks.push_back(_elements.size());
		std::vector<eslocal> ePartition;
		if (blocks.size() == 2) {
			ePartition = getPartition(0, _elements.size(), parts);
		} else {
			ESINFO(ALWAYS) << "NONCONTINUITY" << environment->MPIrank << " (" << blocks.size() - 1 << ").";
			_continuous = false;
			_continuousPartId.clear();
			double averageDomainSize = _elements.size() / (double)parts;
			std::vector<size_t> bPart(blocks.size() - 1);
			size_t bParts = 0;
			for (size_t b = 0; b < blocks.size() - 1; b++) {
				bPart[b] = std::round((blocks[b + 1] - blocks[b]) / averageDomainSize);
				bParts += bPart[b];
			}
			while (bParts < parts) {
				(*std::max_element(bPart.begin(), bPart.end()))++;
				bParts++;
			}
			while (bParts > parts) {
				(*std::max_element(bPart.begin(), bPart.end()))--;
				bParts--;
			}
			for (size_t b = 0; b < bPart.size(); b++) {
				if (bPart[b] == 0) {
					bPart[b]++;
				}
			}
			bParts = 0;
			for (size_t b = 0; b < blocks.size() - 1; b++) {
				std::vector<eslocal> bPartition = getPartition(blocks[b], blocks[b + 1], bPart[b]);
				std::for_each(bPartition.begin(), bPartition.end(), [&] (eslocal &e) { e += bParts; });
				ePartition.insert(ePartition.end(), bPartition.begin(), bPartition.end());
				bParts += bPart[b];
				_continuousPartId.insert(_continuousPartId.end(), bPart[b], b);
			}
			parts = bParts;
		};

		_partPtrs = std::vector<eslocal>(parts + 1, 0);
		for (size_t i = 0; i < _elements.size(); i++) {
			_elements[i]->domains().clear();
			_elements[i]->domains().push_back(ePartition[i]);
			_partPtrs[ePartition[i]]++;
		}

		std::sort(_elements.begin(), _elements.end(), [] (const Element* e1, const Element* e2) { return e1->domains()[0] < e2->domains()[0]; });
		ESTEST(MANDATORY) << "subdomain without element" << (std::any_of(_partPtrs.begin(), _partPtrs.end() - 1, [] (eslocal size) { return size == 0; }) ? TEST_FAILED : TEST_PASSED);
		Esutils::sizesToOffsets(_partPtrs);
	}

	_DOFtoElement.clear();
	_fixPoints.clear();
	_corners.clear();
	mapFacesToClusters();
	mapFacesToDomains();
	mapEdgesToClusters();
	mapEdgesToDomains();
	mapNodesToDomains();
	mapCoordinatesToDomains();
	if (_properties.size()) {
		fillDomainsSettings();
	}
}

void Mesh::partitiate(size_t parts)
{
	ESTEST(MANDATORY) << "Number of domains cannot be " << parts << (parts == 0 ? TEST_FAILED : TEST_PASSED);
	_continuousPartId.clear();
	_continuousPartId.resize(parts, 0);
	if (parts == 1 && this->parts() == 1) {
		_partPtrs = { 0, (eslocal)_elements.size() };
		mapElementsToDomains();
	} else {
		std::vector<eslocal> blocks = continuousReorder(_elements, 0, _elements.size());
		std::vector<eslocal> ePartition;
		if (blocks.size() == 2) {
			ePartition = getPartition(0, _elements.size(), parts);
		} else {
			ESINFO(ALWAYS) << "NONCONTINUITY" << environment->MPIrank << " (" << blocks.size() - 1 << ").";
			_continuous = false;
			_continuousPartId.clear();
			double averageDomainSize = _elements.size() / (double)parts;
			std::vector<size_t> bPart(blocks.size() - 1);
			size_t bParts = 0;
			for (size_t b = 0; b < blocks.size() - 1; b++) {
				bPart[b] = std::round((blocks[b + 1] - blocks[b]) / averageDomainSize);
				bParts += bPart[b];
			}
			while (bParts < parts) {
				(*std::max_element(bPart.begin(), bPart.end()))++;
				bParts++;
			}
			while (bParts > parts) {
				(*std::max_element(bPart.begin(), bPart.end()))--;
				bParts--;
			}
			for (size_t b = 0; b < bPart.size(); b++) {
				if (bPart[b] == 0) {
					bPart[b]++;
				}
			}
			bParts = 0;
			for (size_t b = 0; b < blocks.size() - 1; b++) {
				std::vector<eslocal> bPartition = getPartition(blocks[b], blocks[b + 1], bPart[b]);
				std::for_each(bPartition.begin(), bPartition.end(), [&] (eslocal &e) { e += bParts; });
				ePartition.insert(ePartition.end(), bPartition.begin(), bPartition.end());
				bParts += bPart[b];
				_continuousPartId.insert(_continuousPartId.end(), bPart[b], b);
			}
			parts = bParts;
		};

		_partPtrs = std::vector<eslocal>(parts + 1, 0);
		for (size_t i = 0; i < _elements.size(); i++) {
			_elements[i]->domains().clear();
			_elements[i]->domains().push_back(ePartition[i]);
			_partPtrs[ePartition[i]]++;
		}

		std::sort(_elements.begin(), _elements.end(), [] (const Element* e1, const Element* e2) { return e1->domains()[0] < e2->domains()[0]; });
		ESTEST(MANDATORY) << "subdomain without element" << (std::any_of(_partPtrs.begin(), _partPtrs.end() - 1, [] (eslocal size) { return size == 0; }) ? TEST_FAILED : TEST_PASSED);
		Esutils::sizesToOffsets(_partPtrs);
	}

	_DOFtoElement.clear();
	_fixPoints.clear();
	_corners.clear();
	mapFacesToClusters();
	mapFacesToDomains();
	mapEdgesToClusters();
	mapEdgesToDomains();
	mapNodesToDomains();
	mapCoordinatesToDomains();
	if (_properties.size()) {
		fillDomainsSettings();
	}
}

void APIMesh::partitiate(size_t parts)
{
	ESTEST(MANDATORY) << "Number of domains cannot be " << parts << (parts == 0 ? TEST_FAILED : TEST_PASSED);
	if (_elements.size() / parts < parts / 10.0) {
		parts = _elements.size() / 20;
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "WARNINK: Too small domains. ESPRESO change DOMAINS=" << parts;
	}

	_continuousPartId.clear();
	_continuousPartId.resize(parts, 0);
	if (parts == 1 && this->parts() == 1) {
		_partPtrs = { 0, (eslocal)_elements.size() };
		mapElementsToDomains();
	} else {
		std::vector<eslocal> blocks = continuousReorder(_elements, 0, _elements.size());
		std::vector<eslocal> ePartition;
		if (blocks.size() == 2) {
			ePartition = getPartition(0, _elements.size(), parts);
		} else {
			ESINFO(ALWAYS) << "NONCONTINUITY" << environment->MPIrank << " (" << blocks.size() - 1 << ").";
			_continuous = false;
			_continuousPartId.clear();
			double averageDomainSize = _elements.size() / (double)parts;
			std::vector<size_t> bPart(blocks.size() - 1);
			size_t bParts = 0;
			for (size_t b = 0; b < blocks.size() - 1; b++) {
				bPart[b] = std::round((blocks[b + 1] - blocks[b]) / averageDomainSize);
				bParts += bPart[b];
			}
			while (bParts < parts) {
				(*std::max_element(bPart.begin(), bPart.end()))++;
				bParts++;
			}
			while (bParts > parts) {
				(*std::max_element(bPart.begin(), bPart.end()))--;
				bParts--;
			}
			for (size_t b = 0; b < bPart.size(); b++) {
				if (bPart[b] == 0) {
					bPart[b]++;
				}
			}
			bParts = 0;
			for (size_t b = 0; b < blocks.size() - 1; b++) {
				std::vector<eslocal> bPartition = getPartition(blocks[b], blocks[b + 1], bPart[b]);
				std::for_each(bPartition.begin(), bPartition.end(), [&] (eslocal &e) { e += bParts; });
				ePartition.insert(ePartition.end(), bPartition.begin(), bPartition.end());
				bParts += bPart[b];
				_continuousPartId.insert(_continuousPartId.end(), bPart[b], b);
			}
			parts = bParts;
		};

		_partPtrs = std::vector<eslocal>(parts + 1, 0);
		for (size_t i = 0; i < _elements.size(); i++) {
			_elements[i]->domains().clear();
			_elements[i]->domains().push_back(ePartition[i]);
			_partPtrs[ePartition[i]]++;
		}

		std::sort(_elements.begin(), _elements.end(), [] (const Element* e1, const Element* e2) { return e1->domains()[0] < e2->domains()[0]; });
		ESTEST(MANDATORY) << "subdomain without element" << (std::any_of(_partPtrs.begin(), _partPtrs.end() - 1, [] (eslocal size) { return size == 0; }) ? TEST_FAILED : TEST_PASSED);
		Esutils::sizesToOffsets(_partPtrs);
	}

	_DOFtoElement.clear();
	mapNodesToDomains();
	mapDOFsToDomains();
	mapCoordinatesToDomains();
}

std::vector<eslocal> Mesh::getPartition(const std::vector<Element*> &elements, size_t begin, size_t end, eslocal parts) const
{
	if (parts == 1) {
		return std::vector<eslocal> (end - begin, 0);
	}
	eslocal nCommon;
	std::vector<eslocal> nodeProjection, metisElements, metisNodes;

	METISMeshRepresentation(elements, begin, end, nodeProjection, metisElements, metisNodes, nCommon);

	return computePartition(metisElements, metisNodes, end - begin, nodeProjection.size(), nCommon, parts);
}

std::vector<eslocal> Mesh::getPartition(size_t begin, size_t end, eslocal parts) const
{
	if (parts == 1) {
		return std::vector<eslocal> (end - begin, 0);
	}
	eslocal ncommon;
	std::vector<eslocal> e, n;

	// create array storing pointers to elements' nodes
	e.reserve(end - begin + 1);
	e.push_back(0);
	for (size_t i = begin; i < end; i++) {
		e.push_back(e.back() + _elements[i]->coarseNodes());
	}

	// create array of nodes
	n.reserve(e.back());
	// number of common nodes to be neighbor
	ncommon = 4;
	for (size_t i = begin; i < end; i++) {
		n.insert(n.end(), _elements[i]->indices(), _elements[i]->indices() + _elements[i]->coarseNodes());
		if (ncommon > _elements[i]->nCommon()) {
			ncommon = _elements[i]->nCommon();
		}
	}

	return computePartition(e, n, end - begin, _coordinates->clusterSize(), ncommon, parts);
}

static eslocal computeCenter(std::vector<std::vector<eslocal> > &neighbours)
{
	std::vector<eslocal> ia, ja;

	ia.reserve(neighbours.size() + 1);
	ia.push_back(0);
	for (size_t i = 0; i < neighbours.size(); i++) {
		std::sort(neighbours[i].begin(), neighbours[i].end());
		Esutils::removeDuplicity(neighbours[i]);
		ia.push_back(ia.back() + neighbours[i].size());
		ja.insert(ja.end(), neighbours[i].begin(), neighbours[i].end());
	}

	std::vector<float> a(ja.size(), 1);

	eslocal nSize = neighbours.size();
	std::vector<float> x(nSize, 1. / nSize), y(nSize);

	// Initial vector
	float last_l = nSize, l = 1;
	eslocal incr = 1;

	while (fabs((l - last_l) / l) > 1e-6) {
		mkl_cspblas_scsrsymv("U", &nSize, a.data(), ia.data(), ja.data(), x.data(), y.data());
		last_l = l;
		l = snrm2(&nSize, y.data(), &incr);
		cblas_sscal(nSize, 1 / l, y.data(), incr);
		x.swap(y);
	}

	return cblas_isamax(nSize, x.data(), incr);
}

eslocal Mesh::getCentralNode(const std::vector<Element*> &elements, size_t begin, size_t end, const std::vector<eslocal> &ePartition, eslocal subpart) const
{
	std::vector<eslocal> nMap;
	for (size_t e = begin; e < end; e++) {
		if (ePartition[e - begin] == subpart) {
			for (size_t n = 0; n < elements[e]->nodes(); n++) {
				nMap.push_back(elements[e]->node(n));
			}
		}
	}

	if (!nMap.size()) {
		return -1;
	}

	std::sort(nMap.begin(), nMap.end());
	Esutils::removeDuplicity(nMap);

	auto localIndex = [&] (eslocal index) {
		return std::lower_bound(nMap.begin(), nMap.end(), index) - nMap.begin();
	};

	std::vector<std::vector<eslocal> > neighbours(nMap.size());
	for (size_t e = begin; e < end; e++) {
		if (ePartition[e - begin] == subpart) {
			for (size_t n = 0; n < elements[e]->nodes(); n++) {
				std::vector<eslocal> neigh = elements[e]->getNeighbours(n);
				for (size_t i = 0; i < neigh.size(); i++) {
					if (elements[e]->node(n) < neigh[i]) {
						neighbours[localIndex(elements[e]->node(n))].push_back(localIndex(neigh[i]));
					}
				}
			}
		}
	}

	return nMap[computeCenter(neighbours)];
}

eslocal Mesh::getCentralNode(eslocal begin, eslocal end, const std::vector<eslocal> &ePartition, eslocal part, eslocal subpart) const
{
	// Compute CSR format of symmetric adjacency matrix
	////////////////////////////////////////////////////////////////////////////
	std::vector<std::vector<eslocal> > neighbours(_coordinates->localSize(part));
	for (eslocal i = begin; i < end; i++) {
		if (ePartition[i - begin] == subpart) {
			for (size_t j = 0; j < _elements[i]->nodes(); j++) {
				std::vector<eslocal> neigh = _elements[i]->getNeighbours(j);
				for (size_t k = 0; k < neigh.size(); k++) {
					if (_elements[i]->node(j) < neigh[k]) {
						neighbours[_coordinates->localIndex(_elements[i]->node(j), part)].push_back(_coordinates->localIndex(neigh[k], part));
					}
				}
			}
		}
	}

	return _coordinates->clusterIndex(computeCenter(neighbours), part);
}

Mesh::~Mesh()
{
	for (size_t i = 0; i < _elements.size(); i++) {
		delete _elements[i];
	}

	for (size_t i = 0; i < _faces.size(); i++) {
		delete _faces[i];
	}

	for (size_t i = 0; i < _edges.size(); i++) {
		delete _edges[i];
	}

	for (size_t i = 0; i < _nodes.size(); i++) {
		delete _nodes[i];
	}

	for (size_t i = 0; i < _regions.size(); i++) {
		delete _regions[i];
	}

	for (size_t i = 0; i < _evaluators.size(); i++) {
		delete _evaluators[i];
	}

	for (size_t i = 0; i < _materials.size(); i++) {
		delete _materials[i];
	}
	delete _coordinates;
}

//static bool isOuterFace(
//		std::vector<std::vector<eslocal> > &nodesElements,
//		std::vector<eslocal> &face)
//{
//	std::vector<eslocal> result(nodesElements[face[0]]);
//	std::vector<eslocal>::iterator it = result.end();
//
//	for (size_t i = 1; i < face.size(); i++) {
//		std::vector<eslocal> tmp(result.begin(), it);
//		it = std::set_intersection(tmp.begin(), tmp.end(),
//				nodesElements[face[i]].begin(), nodesElements[face[i]].end(),
//				result.begin());
//		if (it - result.begin() == 1) {
//			return true;
//		}
//	}
//
//	return false;
//}

void Mesh::getSurface(Mesh &surface) const
{
	ESINFO(ERROR) << "Not use getSurface method\n";
//	// vector of faces in all parts
//	std::vector<std::vector<std::vector<eslocal> > > faces(parts());
//	// number of elements in all parts
//	std::vector<size_t> elementsCount(parts(), 0);
//
//	if (parts() < 1) {
//		ESINFO(ERROR) << "Internal error: _partPtrs.size().";
//	}
//
//	#pragma omp parallel for
//	for  (size_t i = 0; i < parts(); i++) {
//		// Compute nodes' adjacent elements
//		std::vector<std::vector<eslocal> > nodesElements(_coordinates->localSize(i));
//		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
//			for (size_t k = 0; k < _elements[j]->nodes(); k++) {
//				nodesElements[_elements[j]->node(k)].push_back(j);
//			}
//		}
//
//		// compute number of elements and fill used nodes
//		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
//			for (size_t k = 0; k < _elements[j]->faces(); k++) {
//				std::vector<eslocal> face(_elements[j]->face(k)->indices(), _elements[j]->face(k)->indices() + _elements[j]->face(k)->nodes());
//				if (isOuterFace(nodesElements, face)) {
//					for (size_t f = 0; f < face.size(); f++) {
//						face[f] = _coordinates->clusterIndex(face[f], i);
//					}
//					faces[i].push_back(face);
//					if (face.size() == 3) {
//						elementsCount[i] += 1;
//					}
//					if (face.size() == 4) {
//						elementsCount[i] += 2;
//					}
//				}
//			}
//		}
//	}
//
//	surface._coordinates = _coordinates;
//
//	size_t count = 0;
//	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
//		count += elementsCount[i];
//	}
//
//	surface._elements.reserve(count);
//	surface._partPtrs.clear();
//	surface._partPtrs.reserve(_partPtrs.size());
//	eslocal params[6] = {0, 0, 0, 0, 0, 0};
//
//	// create surface mesh
//	surface._partPtrs.push_back(0); //(surface._elements.size());
//	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
//		for (size_t j = 0; j < faces[i].size(); j++) {
//			std::vector<eslocal> &face = faces[i][j];
//			if (face.size() == 3) {
//				surface._elements.push_back(new Triangle3(&face[0], params));
//			}
//			// divide square to triangles
//			if (face.size() == 4) {
//				size_t min = 0;
//				for (size_t p = 1; p < 4; p++) {
//					if (_coordinates[face[p]] < _coordinates[face[min]]) {
//						min = p;
//					}
//				}
//				if (min % 2 == 0) {
//					surface._elements.push_back(new Triangle3(&face[0], params));
//					face[1] = face[0];
//					surface._elements.push_back(new Triangle3(&face[1], params));
//				} else {
//					surface._elements.push_back(new Triangle3(&face[1], params));
//					face[2] = face[3];
//					surface._elements.push_back(new Triangle3(&face[0], params));
//				}
//			}
//		}
//		surface._partPtrs.push_back(surface._elements.size());
//	}
//
//	surface.mapCoordinatesToDomains();
//	surface._neighbours = _neighbours;
//	surface._materials = _materials;
//	for (size_t i = 0; i < _evaluators.size(); i++) {
//		surface._evaluators.push_back(_evaluators[i]->copy());
//	}
}

void Mesh::computeVolumeCorners(size_t number, bool onVertices, bool onEdges, bool onFaces)
{
	if (parts() == 1 || (!onVertices && !onEdges && !onFaces)) {
		return;
	}

	if (_corners.size()) {
		return;
	}

	computeFacesSharedByDomains();
	computeEdgesOnBordersOfFacesSharedByDomains();

	computeCornersOnEdges(number, onVertices, onEdges);
	if (onFaces) {
		computeCornersOnFaces(number, onVertices, onEdges, onFaces);
	}
}

void Mesh::computePlaneCorners(size_t number, bool onVertices, bool onEdges)
{
	if (parts() == 1 || (!onVertices && !onEdges)) {
		return;
	}

	if (_corners.size()) {
		return;
	}

	computeEdgesSharedByDomains();
	computeCornersOnEdges(number, onVertices, onEdges);
}

template<typename MergeFunction>
static void uniqueWithMerge(std::vector<Element*> &elements, MergeFunction merge)
{
	if (elements.size() == 0) {
		return;
	}

	auto it = elements.begin();
	auto last = it;
	while (++it != elements.end()) {
		if (!(**last == **it)) {
			*(++last) = *it;
		} else { // Merge two edges
			if (last != it) {
				merge(*last, *it);
				delete *it;
			}
		}
	}
	elements.resize(++last - elements.begin());
}

template<typename MergeFunction>
static std::vector<Element*> mergeElements(size_t threads, std::vector<std::vector<Element*> > &elements, MergeFunction merge)
{
	std::vector<Element*> result;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::sort(elements[t].begin(), elements[t].end(), [] (const Element* e1, const Element* e2) { return *e1 < *e2; });
		uniqueWithMerge(elements[t], merge);
	}

	std::vector<std::vector<Element*> > divided;
	std::vector<std::vector<Element*> > merged;

	divided.swap(elements);
	while (divided.size() > 1) {
		divided.resize(divided.size() + divided.size() % 2); // keep the size even
		merged.resize(divided.size() / 2);

		#pragma omp parallel for
		for (size_t t = 0; t < merged.size(); t++) {
			merged[t].resize(divided[2 * t].size() + divided[2 * t + 1].size());
			std::merge(
					divided[2 * t    ].begin(), divided[2 * t    ].end(),
					divided[2 * t + 1].begin(), divided[2 * t + 1].end(),
					merged[t].begin(), [] (const Element* e1, const Element* e2) { return *e1 < *e2; });
			uniqueWithMerge(merged[t], merge);
		}
		divided.swap(merged);
		merged.clear();
	}

	result.swap(divided[0]);
	return result;
}

const std::vector<Property>& Mesh::propertyGroup(Property property) const
{
	auto it = _propertyGroups.find(property);
	if (it == _propertyGroups.end()) {
		ESINFO(ERROR) << "ESPRESO internal error: request for unknown property group.";
	}
	return it->second;
}

void Mesh::addPropertyGroup(const std::vector<Property> &properties)
{
	for (size_t i = 0; i < properties.size(); i++) {
		_propertyGroups[properties[i]] = properties;
	}
}

void Mesh::loadProperty(
		size_t loadStep,
		const std::map<std::string, std::string> &regions,
		const std::vector<std::string> &parameters,
		const std::vector<Property> &properties,
		ElementType type)
{
	addPropertyGroup(properties);
	auto getValue = [] (const std::vector<std::string> &values, const std::string &parameter) -> std::string {
		for (size_t i = 0; i < values.size(); i++) {
			std::vector<std::string> args = Parser::split(Parser::strip(values[i]), " ");
			if (StringCompare::caseSensitiveEq(args[0], parameter)) {
				std::stringstream ss;
				std::for_each(args.begin() + 1, args.end(), [&] (const std::string &v) { ss << v; });
				return ss.str();
			}
		}
		return "";
	};

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_neighbours.begin(), _neighbours.end(), neighbour) - _neighbours.begin();
	};

	auto distribute = [] (std::vector<Element*> &elements, Property property, Evaluator *evaluator, Region *region) {
		size_t threads = environment->OMP_NUM_THREADS;
		std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				elements[n]->regions().push_back(region);
			}
		}
	};

	auto distributeNodes = [&] (std::vector<eslocal> &nodes, Property property, Evaluator *evaluator, Region *region) {
		size_t threads = environment->OMP_NUM_THREADS;
		std::vector<size_t> distribution = Esutils::getDistribution(threads, nodes.size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				_nodes[nodes[n]]->regions().push_back(region);
			}
		}
	};

	for (auto it = regions.begin(); it != regions.end(); ++it) {
		Region *region = this->region(it->first);
		if (region->settings.size() < loadStep + 1) {
			region->settings.resize(loadStep + 1);
		}
		std::vector<std::string> values = Parser::split(it->second, ",");

		for (size_t p = 0; p < properties.size(); p++) {
			std::string value = properties.size() == 1 ? values[0] : getValue(values, parameters[p]);
			if (!value.size()) {
				continue;
			}

			if (!StringCompare::contains(value, { "x", "y", "z", "TEMPERATURE", "TIME" })) {
				Expression expr(value, {});
				_evaluators.push_back(new ConstEvaluator(expr.evaluate({}), properties[p]));
			} else {
				_evaluators.push_back(new CoordinatesEvaluator(value, *_coordinates, properties[p]));
			}

			std::vector<std::vector<esglobal> > boundaryNodes(_neighbours.size());
			std::vector<std::vector<esglobal> > neighboursNodes(_neighbours.size());
			if (type == ElementType::NODES) {
				ESINFO(OVERVIEW) << "Set " << properties[p] << " to '" << value << "' for LOAD STEP " << loadStep + 1 << " for nodes of region '" << region->name << "'";
			} else if (type == ElementType::FACES) {
				ESINFO(OVERVIEW) << "Set " << properties[p] << " to '" << value << "' for LOAD STEP " << loadStep + 1 << " for faces of region '" << region->name << "'";
			} else {
				ESINFO(OVERVIEW) << "Set " << properties[p] << " to '" << value << "' for LOAD STEP " << loadStep + 1 << " for region '" << region->name << "'";
			}

			if (type == ElementType::NODES && region->elements().size() && region->elements()[0]->nodes() > 1) {
				std::vector<Element*> nodes;
				for (size_t i = 0; i < region->elements().size(); i++) {
					for (size_t n = 0; n < region->elements()[i]->nodes(); n++) {
						nodes.push_back(_nodes[region->elements()[i]->node(n)]);
					}
				}
				std::sort(nodes.begin(), nodes.end());
				Esutils::removeDuplicity(nodes);
				for (size_t n = 0; n < nodes.size(); n++) {
					for (size_t c = 0; c < nodes[n]->clusters().size(); c++) {
						if (nodes[n]->clusters()[c] != environment->MPIrank) {
							boundaryNodes[n2i(nodes[n]->clusters()[c])].push_back(_coordinates->globalIndex(nodes[n]->node(0)));
						}
					}
				}
				distribute(nodes, properties[p], _evaluators.back(), region);
			} else if (type == ElementType::FACES && region->elements().size() && region->elements()[0]->type() != espreso::Element::Type::PLANE) {
				std::vector<std::vector<Element*> > faces(this->parts());
				std::vector<eslocal> nodes;
				for (size_t i = 0; i < region->elements().size(); i++) {
					for (size_t n = 0; n < region->elements()[i]->nodes(); n++) {
						nodes.push_back(region->elements()[i]->node(n));
					}
				}
				std::sort(nodes.begin(), nodes.end());
				Esutils::removeDuplicity(nodes);

				#pragma omp parallel for
				for (size_t d = 0; d < this->parts(); d++) {
					for (eslocal e = this->getPartition()[d]; e < this->getPartition()[d + 1]; e++) {
						Element *face = _elements[e]->addFace(nodes);
						if (face != NULL) {
							faces[d].push_back(face);
						}
					}
				}

				std::vector<Element*> mfaces = mergeElements(this->parts(), faces, [] (Element* e1, Element *e2) {
					e1->parentElements().insert(e1->parentElements().end(), e2->parentElements().begin(), e2->parentElements().end());
					for (size_t e = 0; e < e2->parentElements().size(); e++) {
						for (size_t i = 0; i < e2->parentElements()[e]->faces(); i++) {
							if (e2->parentElements()[0]->face(i) != NULL && *(e1) == *(e2->parentElements()[e]->face(i))) {
								e2->parentElements()[e]->setFace(i, e1);
								break;
							}
						}
					}
				});
				std::vector<std::vector<Element*> > tmp = { _faces, mfaces };

				_faces = mergeElements(2, tmp, [] (Element* e1, Element *e2) {
					e1->parentElements().insert(e1->parentElements().end(), e2->parentElements().begin(), e2->parentElements().end());
					for (size_t e = 0; e < e2->parentElements().size(); e++) {
						for (size_t i = 0; i < e2->parentElements()[e]->faces(); i++) {
							if (e2->parentElements()[0]->face(i) != NULL && *(e1) == *(e2->parentElements()[e]->face(i))) {
								e2->parentElements()[e]->setFace(i, e1);
								break;
							}
						}
					}
				});
				mapFacesToClusters();
				mapFacesToDomains();
				fillParentFacesToNodes();

				region->elements() = mfaces;
				distribute(region->elements(), properties[p], _evaluators.back(), region);
				distributeNodes(nodes, properties[p], _evaluators.back(), region);
			} else {
				std::vector<eslocal> nodes;
				for (size_t i = 0; i < region->elements().size(); i++) {
					for (size_t n = 0; n < region->elements()[i]->nodes(); n++) {
						nodes.push_back(region->elements()[i]->node(n));
					}
				}
				std::sort(nodes.begin(), nodes.end());
				Esutils::removeDuplicity(nodes);

				distribute(region->elements(), properties[p], _evaluators.back(), region);
				distributeNodes(nodes, properties[p], _evaluators.back(), region);
			}

			if (type == ElementType::NODES) {
				if (!Communication::exchangeUnknownSize(boundaryNodes, neighboursNodes, _neighbours)) {
					ESINFO(ERROR) << "problem while synchronization of nodes regions.";
				}
				for (size_t n = 0; n < _neighbours.size(); n++) {
					for (size_t i = 0; i < neighboursNodes[n].size(); i++) {
						_nodes[_coordinates->clusterIndex(neighboursNodes[n][i])]->regions().push_back(region);
					}
				}
			}

			region->settings[loadStep][properties[p]].push_back(_evaluators.back());
		}
	}
}

void Mesh::loadProperty(const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, size_t loadStep)
{
	loadProperty(loadStep, regions, parameters, properties, ElementType::ELEMENTS);
}

void Mesh::loadNodeProperty(const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, size_t loadStep)
{
	loadProperty(loadStep, regions, parameters, properties, ElementType::NODES);
}

void Mesh::loadFaceProperty(const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, size_t loadStep)
{
	loadProperty(loadStep, regions, parameters, properties, ElementType::FACES);
}

void Mesh::loadProperty(const std::map<size_t, std::map<std::string, std::string> > &property, const std::vector<std::string> &parameters, const std::vector<Property> &properties)
{
	for (auto it = property.begin(); it != property.end(); ++it) {
		loadProperty(it->first - 1, it->second, parameters, properties, ElementType::ELEMENTS);
	}
	for (size_t i = 0; i < properties.size(); i++) {
		_propertyGroups[properties[i]] = properties;
	}
}

void Mesh::loadNodeProperty(const std::map<size_t, std::map<std::string, std::string> > &property, const std::vector<std::string> &parameters, const std::vector<Property> &properties)
{
	for (auto it = property.begin(); it != property.end(); ++it) {
		loadProperty(it->first - 1, it->second, parameters, properties, ElementType::NODES);
	}
	for (size_t i = 0; i < properties.size(); i++) {
		_propertyGroups[properties[i]] = properties;
	}
}

void Mesh::loadFaceProperty(const std::map<size_t, std::map<std::string, std::string> > &property, const std::vector<std::string> &parameters, const std::vector<Property> &properties)
{
	for (auto it = property.begin(); it != property.end(); ++it) {
		loadProperty(it->first - 1, it->second, parameters, properties, ElementType::FACES);
	}
	for (size_t i = 0; i < properties.size(); i++) {
		_propertyGroups[properties[i]] = properties;
	}
}

Region* Mesh::region(const std::string &name) const
{
	auto it = std::find_if(_regions.begin(), _regions.end(), [&] (const Region *region) { return region->name.compare(name) == 0; });
	if (it != _regions.end()) {
		return *it;
	}
	ESINFO(GLOBAL_ERROR) << "Unknown region '" << name << "'";
	exit(EXIT_FAILURE);
}

void Mesh::addMonitoredRegion(Region* region) const
{
	for (size_t r = 0; r < _monitoredRegions.size(); r++) {
		if (_monitoredRegions[r] == region) {
			return;
		}
	}
	_monitoredRegions.push_back(region);
}

std::vector<std::vector<Region*> > Mesh::getRegionsWithProperties(const std::vector<Region*> &regions, size_t loadStep, const std::vector<Property> &properties)
{
	std::vector<std::vector<Region*> > result(properties.size());
	for (size_t r = 0; r < regions.size(); r++) {
		Region* region = regions[r];
		for (size_t dof = 0; dof < properties.size(); dof++) {
			if (loadStep < region->settings.size() && region->settings[loadStep].count(properties[dof])) {
				result[dof].push_back(region);
			}
		}
	}
	for (size_t dof = 0; dof < properties.size(); dof++) {
		std::sort(result[dof].begin(), result[dof].end());
		Esutils::removeDuplicity(result[dof]);
	}

	return result;
}

std::vector<std::vector<Region*> > Mesh::getRegionsWithProperties(size_t loadStep, const std::vector<Property> &properties) const
{
	return getRegionsWithProperties(_regions, loadStep, properties);
}

bool Mesh::commonRegion(const std::vector<Region*> &v1, const std::vector<Region*> &v2)
{
	for (size_t i = 0, j = 0; i < v1.size() && j < v2.size(); v1[i] < v2[j] ? i++ : j++) {
		if (v1[i] == v2[j]) {
			return true;
		}
	}
	return false;
}

void Mesh::materialNotFound(const std::string &name)
{
	ESINFO(GLOBAL_ERROR) << "Invalid .ecf file: material " << name << " is not set.";
}

void Mesh::loadMaterial(Region *region, size_t index, const std::string &name, const Configuration &configuration)
{
	#pragma omp parallel for
	for (size_t e = 0; e < region->elements().size(); e++) {
		region->elements()[e]->setParam(Element::MATERIAL, index);
	}
	const Configuration* coordinateSystem = configuration.subconfigurations.find("COORDINATE_SYSTEM")->second;
	_materials.push_back(new Material(*_coordinates, configuration, dynamic_cast<const CoordinateSystem&>(*coordinateSystem)));
	ESINFO(OVERVIEW) << "Set material '" << name << "' for region '" << region->name << "'";
}

void Mesh::checkMaterials()
{
	if (!_materials.size()) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO needs at least one material.";
	}
}

void Mesh::removeDuplicateRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	auto remove = [&] (std::vector<Element*> &elements) {
		std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
				std::sort(elements[i]->regions().begin(), elements[i]->regions().end());
				Esutils::removeDuplicity(elements[i]->regions());
			}
		}
	};

	remove(_elements);
	remove(_faces);
	remove(_edges);
	remove(_nodes);
}

void Mesh::fillDomainsSettings()
{
	_properties.clear();
	_properties.resize(parts());

	for (size_t r = 0; r < _regions.size(); r++) {
		_steps = std::max(_steps, _regions[r]->settings.size());
	}

	auto addProperties = [&] (const std::vector<eslocal> &domains, size_t region) {
		for (size_t d = 0; d < domains.size(); d++) {
			if (_properties[domains[d]].size() < _regions[region]->settings.size()) {
				_properties[domains[d]].resize(_regions[region]->settings.size());
			}
		}
		for (size_t i = 0; i < _regions[region]->settings.size(); i++) {
			for (auto it = _regions[region]->settings[i].begin(); it != _regions[region]->settings[i].end(); ++it) {
				for (size_t d = 0; d < domains.size(); d++) {
					_properties[domains[d]][i].insert(it->first);
				}
			}
		}
	};

	std::vector<eslocal> domains(parts());
	std::iota(domains.begin(), domains.end(), 0);

	addProperties(domains, 0); // settings for ALL_ELEMENTS
	addProperties(domains, 1); // settings for ALL_NODES

	size_t threads = environment->OMP_NUM_THREADS;

	for (size_t r = 2; r < _regions.size(); r++) {
		domains.clear();
		std::vector<size_t> distribution = Esutils::getDistribution(threads, _regions[r]->elements().size());
		std::vector<std::vector<eslocal> > tdomains(threads);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				tdomains[t].insert(tdomains[t].end(), _regions[r]->elements()[e]->domains().begin(), _regions[r]->elements()[e]->domains().end());
			}
			std::sort(tdomains[t].begin(), tdomains[t].end());
			Esutils::removeDuplicity(tdomains[t]);
		}

		for (size_t t = 0; t < threads; t++) {
			domains.insert(domains.end(), tdomains[t].begin(), tdomains[t].end());
		}
		std::sort(domains.begin(), domains.end());
		Esutils::removeDuplicity(domains);
		addProperties(domains, r);
	}
}

bool Mesh::hasProperty(size_t domain, Property property, size_t loadStep) const
{
	return loadStep < _properties[domain].size() && _properties[domain][loadStep].count(property);
}

bool Mesh::hasProperty(Property property, size_t loadStep) const
{
	for (size_t r = 0; r < _regions.size(); r++) {
		if (loadStep < _regions[r]->settings.size() && _regions[r]->settings[loadStep].find(property) != _regions[r]->settings[loadStep].end()) {
			return true;
		}
	}
	return false;
}

bool Mesh::isPropertyTimeDependent(Property property, size_t loadStep) const
{
	for (size_t r = 0; r < _regions.size(); r++) {
		if (loadStep < _regions[r]->settings.size()) {
			auto it = _regions[r]->settings[loadStep].find(property);
			if (it != _regions[r]->settings[loadStep].end()) {
				return std::any_of(it->second.begin(), it->second.end(), [] (const Evaluator *e) { return e->isTimeDependent(); });
			}
		}
	}
	return false;
}

bool Mesh::isPropertyTemperatureDependent(Property property, size_t loadStep) const
{
	for (size_t r = 0; r < _regions.size(); r++) {
		if (loadStep < _regions[r]->settings.size()) {
			auto it = _regions[r]->settings[loadStep].find(property);
			if (it != _regions[r]->settings[loadStep].end()) {
				return std::any_of(it->second.begin(), it->second.end(), [] (const Evaluator *e) { return e->isTemperatureDependent(); });
			}
		}
	}
	return false;
}

bool Mesh::isAnyPropertyTimeDependent(const std::vector<Property> &properties, size_t loadStep) const
{
	return std::any_of(properties.begin(), properties.end(), [&] (Property p) { return isPropertyTimeDependent(p, loadStep); });
}

bool Mesh::isAnyPropertyTemperatureDependent(const std::vector<Property> &properties, size_t loadStep) const
{
	return std::any_of(properties.begin(), properties.end(), [&] (Property p) { return isPropertyTemperatureDependent(p, loadStep); });
}

void Mesh::fillEdgesFromElements(std::function<bool(const std::vector<Element*> &nodes, const Element* edge)> filter)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _elements.size());

	std::vector<std::vector<Element*> > edges(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			size_t i = _elements[e]->fillEdges();

			for (; i < _elements[e]->edges(); i++) {
				Element* edge = _elements[e]->edge(i);
				if (edge->regions().size() || filter(_nodes, edge)) {
					edges[t].push_back(edge);
				} else {
					_elements[e]->setEdge(i, NULL);
					delete edge;
				}
			}
		}
	}

	std::vector<Element*> created = mergeElements(threads, edges, [] (Element* e1, Element *e2) {
		e1->parentElements().insert(e1->parentElements().end(), e2->parentElements().begin(), e2->parentElements().end());
		for (size_t e = 0; e < e2->parentElements().size(); e++) {
			for (size_t i = 0; i < e2->parentElements()[e]->edges(); i++) {
				if (e2->parentElements()[0]->edge(i) != NULL && *(e1) == *(e2->parentElements()[e]->edge(i))) {
					e2->parentElements()[e]->setEdge(i, e1);
					break;
				}
			}
		}
	});
	_edges.insert(_edges.end(), created.begin(), created.end());

	mapEdgesToClusters();
	mapEdgesToDomains();
	fillParentEdgesToNodes();
}

void Mesh::fillFacesFromElements(std::function<bool(const std::vector<Element*> &nodes, const Element* face)> filter)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _elements.size());

	std::vector<std::vector<Element*> > faces(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			size_t f = _elements[e]->fillFaces();

			for (; f < _elements[e]->faces(); f++) {
				Element* face = _elements[e]->face(f);
				if (face->regions().size() || filter(_nodes, face)) {
					faces[t].push_back(face);
				} else {
					_elements[e]->setFace(f, NULL);
					delete face;
				}
			}
		}
	}

	std::vector<Element*> created = mergeElements(threads, faces, [] (Element* e1, Element *e2) {
		e1->parentElements().push_back(e2->parentElements().back()); // Face can be only between two elements
		for (size_t i = 0; i < e2->parentElements()[0]->faces(); i++) {
			if (e2->parentElements()[0]->face(i) != NULL && *e1 == *e2->parentElements()[0]->face(i)) {
				e2->parentElements()[0]->setFace(i, e1);
				break;
			}
		}
	});
	_faces.insert(_faces.end(), created.begin(), created.end());

	mapFacesToClusters();
	mapFacesToDomains();
	fillParentFacesToNodes();
}

void Mesh::fillNodesFromCoordinates()
{
	_nodes.reserve(_coordinates->clusterSize());
	for (size_t i = 0; i < _coordinates->clusterSize(); i++) {
		_nodes.push_back(new Node(i));
	}
}

void Mesh::computeElementsFromFaces()
{
	ESINFO(GLOBAL_ERROR) << "Implement computeElementFromFaces";
}

void Mesh::fillParentElementsToNodes()
{
	#pragma omp parallel for
	for  (size_t i = 0; i < _nodes.size(); i++) {
		_nodes[i]->parentElements().clear();
	}

	for (size_t e = 0; e < _elements.size(); e++) {
		for (size_t n = 0; n < _elements[e]->nodes(); n++) {
			_nodes[_elements[e]->node(n)]->parentElements().push_back(_elements[e]);
		}
	}

	#pragma omp parallel for
	for  (size_t i = 0; i < _nodes.size(); i++) {
		std::sort(_nodes[i]->parentElements().begin(), _nodes[i]->parentElements().end());
	}
}

void APIMesh::fillParentElementsToDOFs(const std::vector<std::vector<eslocal> > &eDOFs)
{
	ESTEST(MANDATORY) << "Invalid number of recognized elements in API." << (eDOFs.size() != _elements.size() ? TEST_FAILED : TEST_PASSED);

	#pragma omp parallel for
	for  (size_t i = 0; i < _DOFs.size(); i++) {
		_DOFs[i]->parentElements().clear();
	}

	for (size_t e = 0; e < eDOFs.size(); e++) {
		for (size_t d = 0; d < eDOFs[e].size(); d++) {
			_DOFs[eDOFs[e][d]]->parentElements().push_back(_elements[e]);
		}
	}

	#pragma omp parallel for
	for  (size_t i = 0; i < _DOFs.size(); i++) {
		std::sort(_DOFs[i]->parentElements().begin(), _DOFs[i]->parentElements().end());
	}
}

void Mesh::fillParentFacesToNodes()
{
	#pragma omp parallel for
	for  (size_t i = 0; i < _nodes.size(); i++) {
		_nodes[i]->parentFaces().clear();
	}

	for (size_t f = 0; f < _faces.size(); f++) {
		for (size_t n = 0; n < _faces[f]->nodes(); n++) {
			_nodes[_faces[f]->node(n)]->parentFaces().push_back(_faces[f]);
		}
	}

	#pragma omp parallel for
	for  (size_t i = 0; i < _nodes.size(); i++) {
		std::sort(_nodes[i]->parentFaces().begin(), _nodes[i]->parentFaces().end());
	}
}

void Mesh::fillParentEdgesToNodes()
{
	#pragma omp parallel for
	for  (size_t i = 0; i < _nodes.size(); i++) {
		_nodes[i]->parentEdges().clear();
	}

	for (size_t e = 0; e < _edges.size(); e++) {
		for (size_t n = 0; n < _edges[e]->nodes(); n++) {
			_nodes[_edges[e]->node(n)]->parentEdges().push_back(_edges[e]);
		}
	}

	#pragma omp parallel for
	for  (size_t i = 0; i < _nodes.size(); i++) {
		std::sort(_nodes[i]->parentEdges().begin(), _nodes[i]->parentEdges().end());
	}
}

void Mesh::fillEdgesFromFaces(std::function<bool(const std::vector<Element*> &faces, const Element* edge)> filter)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());

	std::vector<std::vector<Element*> > edges(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			_faces[e]->fillEdges();

			for (size_t i = 0; i < _faces[e]->edges(); i++) {
				Element* edge = _faces[e]->edge(i);
				if (edge->regions().size() || filter(_nodes, edge)) {
					edges[t].push_back(edge);
				} else {
					_faces[e]->setEdge(i, NULL);
					delete edge;
				}
			}
		}
	}

	_edges = mergeElements(threads, edges, [] (Element* e1, Element *e2) {
		e1->parentElements().insert(e1->parentElements().end(), e2->parentElements().begin(), e2->parentElements().end());
		for (size_t e = 0; e < e2->parentElements().size(); e++) {
			for (size_t i = 0; i < e2->parentElements()[e]->edges(); i++) {
				if (e2->parentElements()[e]->edge(i) != NULL && *(e1) == *(e2->parentElements()[e]->edge(i))) {
					e2->parentElements()[e]->setEdge(i, e1);
					break;
				}
			}
		}
	});

	mapEdgesToClusters();
	mapEdgesToDomains();
	fillParentEdgesToNodes();
}

static Element* parentElement(const std::vector<Element*> &nodes, const Element *e)
{
	std::vector<Element*> intersection(nodes[e->node(e->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
	auto it = intersection.end();

	for (size_t n = e->nodes() - 2; it - intersection.begin() > 1 &&  n < e->nodes(); n--) {
		std::vector<Element*> tmp(intersection.begin(), it);
		it = std::set_intersection(tmp.begin(), tmp.end(),
				nodes[e->node(n)]->parentElements().begin(), nodes[e->node(n)]->parentElements().end(),
				intersection.begin());
	}

	intersection.resize(it - intersection.begin());
	return intersection[0];
}

void Mesh::fillEdgesParents()
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _edges.size());
	std::vector<std::vector<Element*> > edges(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			if (!_edges[e]->parentElements().size()) {
				_edges[e]->addParent(parentElement(_nodes, _edges[e]));
				edges[t].push_back(_edges[e]);
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		for (size_t e = 0; e < edges[t].size(); e++) {
			edges[t][e]->parentElements()[0]->addEdge(edges[t][e]);
		}
	}
}

void Mesh::fillFacesParents()
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());
	std::vector<std::vector<Element*> > faces(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			if (!_faces[e]->parentElements().size()) {
				_faces[e]->addParent(parentElement(_nodes, _faces[e]));
				faces[t].push_back(_faces[e]);
			}
		}
	}

	for (size_t t = 0; t < threads; t++) {
		for (size_t f = 0; f < faces[t].size(); f++) {
			faces[t][f]->parentElements()[0]->addFace(faces[t][f]);
		}
	}
}

void Mesh::computeFacesOfAllElements()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element* face) { return true; });
}


void Mesh::computeFacesOnDomainsSurface()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element *face) {
		std::vector<Element*> intersection(nodes[face->node(face->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = face->nodes() - 2; it - intersection.begin() > 1 &&  n < face->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[face->node(n)]->parentElements().begin(), nodes[face->node(n)]->parentElements().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1) {
			return true;
		}
		return intersection[0]->domains() != intersection[1]->domains();
	});
}

void Mesh::computeFacesSharedByDomains()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element *face) {
		std::vector<Element*> intersection(nodes[face->node(face->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = face->nodes() - 2; it - intersection.begin() > 1 &&  n < face->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[face->node(n)]->parentElements().begin(), nodes[face->node(n)]->parentElements().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1) {
			return false;
		}
		return intersection[0]->domains() != intersection[1]->domains();
	});
}

void APIMesh::computeFacesSharedByDomains()
{
	fillFacesFromElements([] (const std::vector<Element*> &nodes, const Element* face) { return true; });

	#pragma omp parallel for
	for  (size_t f = 0; f < _faces.size(); f++) {
		std::vector<eslocal> d0 = _faces[f]->parentElements()[0]->DOFsIndices();
		std::vector<eslocal> d1 = _faces[f]->parentElements()[1]->DOFsIndices();
		std::sort(d0.begin(), d0.end());
		std::sort(d1.begin(), d1.end());
		std::vector<eslocal> intersection(_faces[f]->parentElements()[0]->DOFsIndices().size());
		auto it = std::set_intersection(d0.begin(), d0.end(), d1.begin(), d1.end(), intersection.begin());
		_faces[f]->DOFsIndices().insert(_faces[f]->DOFsIndices().end(), intersection.begin(), it);
	}
}

void Mesh::clearFacesWithoutSettings()
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t f = distribution[t]; f < distribution[t + 1]; f++) {
			if (!_faces[f]->regions().size()) {
				for (size_t e = 0; e < _faces[f]->parentElements().size(); e++) {
					for (size_t i = 0; i < _faces[f]->parentElements()[e]->faces(); i++) {
						if (_faces[f]->parentElements()[e]->face(i) != NULL && *_faces[f] == *_faces[f]->parentElements()[e]->face(i)) {
							_faces[f]->parentElements()[e]->setFace(i, NULL);
						}
					}
				}
				delete _faces[f];
				_faces[f] = NULL;
			}
		}
	}

	size_t it = 0;
	for (size_t i = 0; i < _faces.size(); i++) {
		if (_faces[i] != NULL) {
			_faces[it++] = _faces[i];
		}
	}
	_faces.resize(it);
}

void Mesh::computeEdgesOfAllElements()
{
	fillEdgesFromElements([] (const std::vector<Element*> &nodes, const Element* edge) { return true; });
}

void Mesh::computeEdgesSharedByDomains()
{
	fillEdgesFromElements([] (const std::vector<Element*> &nodes, const Element *edge) {
		std::vector<Element*> intersection(nodes[edge->node(edge->nodes() - 1)]->parentElements()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = edge->nodes() - 2; it - intersection.begin() > 1 &&  n < edge->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[edge->node(n)]->parentElements().begin(), nodes[edge->node(n)]->parentElements().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1) {
			return false;
		}
		return intersection[0]->domains() != intersection[1]->domains();
	});
}

void Mesh::computeEdgesOnBordersOfFacesSharedByDomains()
{
	fillEdgesFromFaces([] (const std::vector<Element*> &nodes, const Element* edge) {
		std::vector<Element*> intersection(nodes[edge->node(edge->nodes() - 1)]->parentFaces()); // it is better to start from end (from mid points)
		auto it = intersection.end();

		for (size_t n = edge->nodes() - 2; it - intersection.begin() > 1 &&  n < edge->nodes(); n--) {
			std::vector<Element*> tmp(intersection.begin(), it);
			it = std::set_intersection(tmp.begin(), tmp.end(),
					nodes[edge->node(n)]->parentFaces().begin(), nodes[edge->node(n)]->parentFaces().end(),
					intersection.begin());
		}

		intersection.resize(it - intersection.begin());
		if (intersection.size() == 1 && intersection[0]->domains().size() == 2) {
			return true;
		}

		for (size_t i = 1; i < intersection.size(); i++) {
			if (intersection[i]->domains() != intersection[i - 1]->domains()) {
				return true;
			}
		}
		return false;
	});
}

void Mesh::clearEdgesWithoutSettings()
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _edges.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t f = distribution[t]; f < distribution[t + 1]; f++) {
			if (!_edges[f]->regions().size()) {
				for (size_t e = 0; e < _edges[f]->parentElements().size(); e++) {
					for (size_t i = 0; i < _edges[f]->parentElements()[e]->edges(); i++) {
						if (_edges[f]->parentElements()[e]->edge(i) != NULL && *_edges[f] == *_edges[f]->parentElements()[e]->edge(i)) {
							_edges[f]->parentElements()[e]->setEdge(i, NULL);
						}
					}
				}
				delete _edges[f];
				_edges[f] = NULL;
			}
		}
	}

	size_t it = 0;
	for (size_t i = 0; i < _edges.size(); i++) {
		if (_edges[i] != NULL) {
			_edges[it++] = _edges[i];
		}
	}
	_edges.resize(it);
}

void Mesh::computeCornersOnEdges(size_t number, bool onVertices, bool onEdges)
{
	if (_edges.size() == 0) {
		ESINFO(ERROR) << "There are no edges for computation of corners.";
	}

	std::vector<Element*> edges;
	if (onEdges) {
		edges.reserve(_edges.size());
	}
	for (size_t e = 0; e < _edges.size(); e++) {
		if (_edges[e]->domains().size()) {
			if (onVertices) {
				if (_nodes[_edges[e]->node(0)]->domains().size() > _nodes[_edges[e]->node(_edges[e]->coarseNodes() - 1)]->domains().size()) {
					_corners.push_back(_nodes[_edges[e]->node(0)]);
				}
				if (_nodes[_edges[e]->node(0)]->domains().size() < _nodes[_edges[e]->node(_edges[e]->coarseNodes() - 1)]->domains().size()) {
					_corners.push_back(_nodes[_edges[e]->node(_edges[e]->coarseNodes() - 1)]);
				}
			}
			if (onEdges) {
				edges.push_back(_edges[e]);
			}
		}
	}

	if (!edges.size()) {
		return;
	}

	std::sort(edges.begin(), edges.end(), [] (Element* e1, Element* e2) { return e1->domains() < e2->domains(); });

	std::vector<size_t> partPtrs;
	partPtrs.push_back(0);
	for (size_t i = 1; i < edges.size(); i++) {
		if (edges[i]->domains() != edges[i - 1]->domains()) {
			partPtrs.push_back(i);
		}
	}
	partPtrs.push_back(edges.size());

	std::vector<std::vector<Element*> > corners(partPtrs.size() - 1);
	#pragma omp parallel for
	for  (size_t p = 0; p < partPtrs.size() - 1; p++) {
		std::vector<eslocal> subPartPtrs = continuousReorder(edges, partPtrs[p], partPtrs[p + 1]);
		for (size_t sp = 0; sp < subPartPtrs.size() - 1; sp++) {
			std::vector<eslocal> nodes;
			for (eslocal e = subPartPtrs[sp]; e < subPartPtrs[sp + 1]; e++) {
				nodes.push_back(edges[e]->node(0));
				nodes.push_back(edges[e]->node(edges[e]->coarseNodes() - 1));
			}
			std::sort(nodes.begin(), nodes.end());
			bool cycle = nodes.size() % 2 == 0;
			for (size_t i = 0; cycle && i < nodes.size(); i += 2) {
				cycle = nodes[i + 1] == nodes[i];
			}
			if (cycle && nodes.size() < 12) {
				Esutils::removeDuplicity(nodes);
				for (size_t i = 0; i < nodes.size() && i < 4; i++) {
					corners[p].push_back(_nodes[nodes[i]]);
				}
				continue;
			}
			if (!cycle && !number) {
				continue;
			}
			std::vector<eslocal> ePartition = getPartition(edges, subPartPtrs[sp], subPartPtrs[sp + 1], cycle ? 4 : number);
			for (size_t c = 0; c < (cycle ? 4 : number); c++) {
				eslocal center = getCentralNode(edges, subPartPtrs[sp], subPartPtrs[sp + 1], ePartition, c);
				if (center > -1) {
					corners[p].push_back(_nodes[center]);
				}
			}
		}
	}

	for (size_t p = 0; p < corners.size(); p++) {
		_corners.insert(_corners.end(), corners[p].begin(), corners[p].end());
	}
	std::sort(_corners.begin(), _corners.end());
	Esutils::removeDuplicity(_corners);
}

void Mesh::computeCornersOnFaces(size_t number, bool onVertices, bool onEdges, bool onFaces)
{
	ESINFO(GLOBAL_ERROR) << "Corners in faces are not implemented.";
}

static void setCluster(Element* &element, std::vector<Element*> &nodes)
{
	std::vector<eslocal> intersection = nodes[element->node(0)]->clusters();
	for (size_t i = 1; i < element->nodes(); i++) {
		auto tmp(intersection);
		auto it = std::set_intersection(
				nodes[element->node(i)]->clusters().begin(), nodes[element->node(i)]->clusters().end(),
				tmp.begin(), tmp.end(), intersection.begin());
		intersection.resize(it - intersection.begin());
	}
	element->clusters() = intersection;
}

void Mesh::mapFacesToClusters()
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _faces.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t f = distribution[t]; f < distribution[t + 1]; f++) {
			if (_faces[f]->parentElements().size() == 1) { // Only faces with one element can have more clusters
				if (std::all_of(_faces[f]->indices(), _faces[f]->indices() + _faces[f]->coarseNodes(), [&] (eslocal i) { return _nodes[i]->clusters().size() > 1; })) {
					setCluster(_faces[f], _nodes);
					continue;
				}
			}
			_faces[f]->clusters() = { environment->MPIrank };
		}
	}
}

void Mesh::mapEdgesToClusters()
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _edges.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
			if (std::all_of(_edges[e]->indices(), _edges[e]->indices() + _edges[e]->coarseNodes(), [&] (eslocal i) { return _nodes[i]->clusters().size() > 1; })) {
				setCluster(_edges[e], _nodes);
			} else {
				_edges[e]->clusters() = { environment->MPIrank };
			}
		}
	}
}

static void assignDomains(std::vector<Element*> &elements)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			elements[i]->domains().clear();
			for (size_t e = 0; e < elements[i]->parentElements().size(); e++) {
				elements[i]->domains().insert(elements[i]->domains().end(), elements[i]->parentElements()[e]->domains().begin(), elements[i]->parentElements()[e]->domains().end());
			}
			std::sort(elements[i]->domains().begin(), elements[i]->domains().end());
			Esutils::removeDuplicity(elements[i]->domains());
		}
	}
}

void Mesh::mapElementsToDomains()
{
	#pragma omp parallel for
	for  (size_t p = 0; p < parts(); p++) {
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			_elements[e]->domains().clear();
			_elements[e]->domains().push_back(p);
		}
	}
}

void Mesh::mapFacesToDomains()
{
	assignDomains(_faces);
}

void Mesh::mapEdgesToDomains()
{
	assignDomains(_edges);
}

void Mesh::mapNodesToDomains()
{
	assignDomains(_nodes);
}

void APIMesh::mapDOFsToDomains()
{
	assignDomains(_DOFs);
}

static void setDOFsIndices(
		std::vector<Element*> &elements,
		std::vector<std::vector<Element*> > &DOFtoElement,
		size_t parts,
		size_t DOFsSize,
		const std::vector<size_t> &DOFs,
		const std::vector<size_t> &offsets,
		const std::vector<std::vector<std::vector<size_t> > > &threadsOffsets)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	std::vector<std::vector<std::vector<Element*> > > tDOFtoElement(threads, std::vector<std::vector<Element*> >(parts));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<size_t> counters(parts);
		for (size_t p = 0; p < parts; p++) {
			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				counters[p] += threadsOffsets[p][dof][t];
			}
			counters[p] += offsets[p];
		}

		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			for (size_t d = 0; d < elements[i]->domains().size(); d++) {

				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					if (elements[i]->DOFsIndices()[d * DOFsSize + DOFs[dof]] == __HAS_DOF__) {

						elements[i]->DOFsIndices()[d * DOFsSize + DOFs[dof]] = counters[elements[i]->domains()[d]]++;
						tDOFtoElement[t][elements[i]->domains()[d]].push_back(elements[i]);

					}
				}
			}
		}
	}

	DOFtoElement.resize(parts);
	#pragma omp parallel for
	for (size_t p = 0; p < parts; p++) {
		for (size_t t = 0; t < threads; t++) {
			DOFtoElement[p].insert(DOFtoElement[p].end(), tDOFtoElement[t][p].begin(), tDOFtoElement[t][p].end());
		}
	}
}

static void fillDOFsOffsets(std::vector<Property> &allDOFsOffsets, const std::vector<Property> &DOFs, std::vector<size_t> &DOFsOffsets)
{
	DOFsOffsets.resize(DOFs.size());
	for (size_t i = 0; i < DOFs.size(); i++) {
		DOFsOffsets[i] = allDOFsOffsets.size();
		allDOFsOffsets.push_back(DOFs[i]);
	}
}

void Mesh::repartitiate(size_t parts, std::vector<size_t> &domainDOFCount, std::vector<std::vector<eslocal> > &previousDOFMap, std::vector<std::vector<eslocal> > &previousDomainMap)
{
	previousDOFMap.clear();
	previousDomainMap.clear();
	for (size_t n = 0; n < _nodes.size(); n++) {
		previousDOFMap.push_back(_nodes[n]->DOFsIndices());
		previousDomainMap.push_back(_nodes[n]->domains());
		_nodes[n]->DOFsIndices().clear();
	}

	partitiate(parts);

	// TODO: DOFs in edges, faces, elements
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _nodes.size());

	size_t DOFSize = _nodesDOFsOffsets.size();
	std::vector<size_t> offsets(parts);
	std::vector<size_t> DOFs(DOFSize);
	std::iota(DOFs.begin(), DOFs.end(), 0);

	// domains x DOFs x (threads + 1)
	std::vector<std::vector<std::vector<size_t> > > threadsOffsets(parts, std::vector<std::vector<size_t> >(DOFSize, std::vector<size_t>(threads + 1)));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<size_t> > threadOffset(parts, std::vector<size_t>(DOFSize, 0));
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			_nodes[i]->DOFsIndices().reserve(DOFSize * _nodes[i]->domains().size());

			_nodes[i]->DOFsIndices().insert(_nodes[i]->DOFsIndices().begin(), _nodes[i]->domains().size() * DOFSize, __WITHOUT_DOF__);
			for (size_t d = 0; d < _nodes[i]->domains().size(); d++) {
				for (size_t dof = 0; dof < DOFSize; dof++) {
					if (previousDOFMap[i][dof] != __WITHOUT_DOF__) {
						_nodes[i]->DOFsIndices()[d * DOFSize + dof] = __HAS_DOF__;
						threadOffset[_nodes[i]->domains()[d]][dof]++;
					}
				}
			}
		}

		for (size_t p = 0; p < parts; p++) {
			for (size_t d = 0; d < DOFSize; d++) {
				threadsOffsets[p][d][t] = threadOffset[p][d];
			}
		}
	}

	std::vector<size_t> sizes(parts);
	for (size_t p = 0; p < parts; p++) {
		for (size_t d = 0; d < DOFSize; d++) {
			sizes[p] += Esutils::sizesToOffsets(threadsOffsets[p][d]);
		}
	}
	domainDOFCount = sizes;

	_DOFtoElement.clear();
	setDOFsIndices(_nodes, _DOFtoElement, parts, DOFSize, DOFs, offsets, threadsOffsets);
}

std::vector<size_t> Mesh::assignVariousDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs, std::vector<size_t> &DOFsOffsets)
{
	auto findDOF = [&] (const std::vector<Property> &DOFs, size_t &added, std::vector<bool> &addDOF) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (!addDOF[dof] && std::find(DOFs.begin(), DOFs.end(), DOFs[dof]) != DOFs.end()) {
				added++;
				addDOF[dof] = true;
			}
		}
	};

	auto fillDOFs = [&] (Element *node, eslocal i, eslocal domain, std::vector<bool> &addDOF) {
		for (size_t e = 0, added = 0; e < node->parentElements().size() && added < DOFs.size(); e++) {
			const Element *el = node->parentElements()[e];
			if (el->domains()[0] == domain) {
				if (std::find(el->indices(), el->indices() + el->coarseNodes(), i) == el->indices() + el->coarseNodes()) {
					findDOF(el->midPointDOFs(), added, addDOF);
				} else {
					findDOF(el->pointDOFs(), added, addDOF);
				}
			}
		}
	};

	size_t prevDOFsSize = _nodesDOFsOffsets.size();
	fillDOFsOffsets(_nodesDOFsOffsets, DOFs, DOFsOffsets);

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _nodes.size());

	// domains x DOFs x (threads + 1)
	std::vector<std::vector<std::vector<size_t> > > threadsOffsets(parts(), std::vector<std::vector<size_t> >(DOFs.size(), std::vector<size_t>(threads + 1)));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<size_t> > threadOffset(parts(), std::vector<size_t>(DOFs.size(), 0));
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			_nodes[i]->DOFsIndices().reserve(_nodesDOFsOffsets.size() * _nodes[i]->domains().size());

			for (size_t d = 0; d < _nodes[i]->domains().size(); d++) {
				_nodes[i]->DOFsIndices().insert(_nodes[i]->DOFsIndices().begin() + (d + 1) * prevDOFsSize + d * (prevDOFsSize + DOFs.size()), DOFs.size(), __WITHOUT_DOF__);
				std::vector<bool> addDOF(DOFs.size(), false);

				fillDOFs(_nodes[i], i, _nodes[i]->domains()[d], addDOF);

				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					if (addDOF[dof]) {
						_nodes[i]->DOFsIndices()[d * DOFs.size() + dof] = __HAS_DOF__;
						threadOffset[_nodes[i]->domains()[d]][dof]++;
					}
				}
			}
		}

		for (size_t p = 0; p < parts(); p++) {
			for (size_t d = 0; d < DOFs.size(); d++) {
				threadsOffsets[p][d][t] = threadOffset[p][d];
			}
		}
	}

	std::vector<size_t> sizes(offsets);
	for (size_t p = 0; p < parts(); p++) {
		for (size_t d = 0; d < DOFs.size(); d++) {
			sizes[p] += Esutils::sizesToOffsets(threadsOffsets[p][d]);
		}
	}

	setDOFsIndices(_nodes, _DOFtoElement, parts(), _nodesDOFsOffsets.size(), DOFsOffsets, offsets, threadsOffsets);

	return sizes;
}


static std::vector<size_t> fillUniformDOFs(
		std::vector<Element*> &elements,
		std::vector<std::vector<Element*> > &DOFtoElement,
		size_t parts,
		size_t prevDOFsSize,
		const std::vector<size_t> &DOFs,
		const std::vector<size_t> &offsets)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	// domains x DOF x (threads + 1)
	std::vector<std::vector<std::vector<size_t> > > threadsOffsets(parts, std::vector<std::vector<size_t> >(DOFs.size(), std::vector<size_t>(threads + 1)));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<size_t> threadOffset(parts, 0);
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			elements[i]->DOFsIndices().reserve(elements[i]->DOFsIndices().size() + DOFs.size() * elements[i]->domains().size());

			for (size_t d = 0; d < elements[i]->domains().size(); d++) {
				elements[i]->DOFsIndices().insert(elements[i]->DOFsIndices().begin() + (d + 1) * prevDOFsSize + d * DOFs.size(), DOFs.size(), __HAS_DOF__);
				threadOffset[elements[i]->domains()[d]]++;
			}
		}

		for (size_t p = 0; p < parts; p++) {
			for (size_t d = 0; d < DOFs.size(); d++) {
				threadsOffsets[p][d][t] = threadOffset[p];
			}
		}
	}

	std::vector<size_t> sizes(offsets);
	for (size_t p = 0; p < parts; p++) {
		for (size_t d = 0; d < DOFs.size(); d++) {
			sizes[p] += Esutils::sizesToOffsets(threadsOffsets[p][d]);
		}
	}

	setDOFsIndices(elements, DOFtoElement, parts, prevDOFsSize + DOFs.size(), DOFs, offsets, threadsOffsets);

	return sizes;
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs, std::vector<size_t> &DOFsOffsets)
{
	fillDOFsOffsets(_nodesDOFsOffsets, DOFs, DOFsOffsets);
	return fillUniformDOFs(_nodes, _DOFtoElement, parts(), _nodesDOFsOffsets.size() - DOFs.size(), DOFsOffsets, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToEdges(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs, std::vector<size_t> &DOFsOffsets)
{
	fillDOFsOffsets(_edgesDOFsOffsets, DOFs, DOFsOffsets);
	return fillUniformDOFs(_edges, _DOFtoElement, parts(), 0, DOFsOffsets, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToFaces(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs, std::vector<size_t> &DOFsOffsets)
{
	fillDOFsOffsets(_facesDOFsOffsets, DOFs, DOFsOffsets);
	return fillUniformDOFs(_faces, _DOFtoElement, parts(), 0, DOFsOffsets, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToElements(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs, std::vector<size_t> &DOFsOffsets)
{
	fillDOFsOffsets(_elementsDOFsOffsets, DOFs, DOFsOffsets);
	return fillUniformDOFs(_elements, _DOFtoElement, parts(), 0, DOFsOffsets, offsets);
}

std::vector<size_t> APIMesh::distributeDOFsToDomains(const std::vector<size_t> &offsets)
{
	return fillUniformDOFs(_DOFs, _DOFtoElement, parts(), 0, { 0 }, offsets);
}

static void computeDOFsCounters(std::vector<Element*> &elements, const std::vector<Property> &DOFs, std::vector<int> neighbours, const std::vector<esglobal> &l2g, const std::vector<G2L> &g2l)
{
	neighbours.push_back(environment->MPIrank);
	std::sort(neighbours.begin(), neighbours.end());

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());

	// threads x neighbour x data
	std::vector<std::vector<std::vector<esglobal> > > sBuffer(threads, std::vector<std::vector<esglobal> >(neighbours.size()));
	// neighbour x data
	std::vector<std::vector<esglobal> > rBuffer(neighbours.size());

	size_t prevDOFsSize = 0;
	if (elements.size()) {
		prevDOFsSize = elements[0]->DOFsDomainsCounters().size() / elements[0]->clusters().size();
	}

	// Compute send buffers
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (elements[e]->nodes() == 1 && elements[e]->parentElements().size() == 0) {
				// mesh generator can generate dangling nodes -> skip them
				elements[e]->clusters().clear();
				continue;
			}

			eslocal cluster = 0;
			elements[e]->DOFsDomainsCounters().reserve(elements[e]->DOFsDomainsCounters().size() + DOFs.size() * elements[e]->clusters().size());
			elements[e]->numberOfGlobalDomains(elements[e]->domains().size());
			for (size_t c = 0; c < elements[e]->clusters().size(); c++) {
				elements[e]->DOFsDomainsCounters().insert(elements[e]->DOFsDomainsCounters().begin() + (c + 1) * prevDOFsSize + c * DOFs.size(), DOFs.size(), -1);
				if (elements[e]->clusters()[c] == environment->MPIrank) {
					cluster = c;
					for (size_t i = prevDOFsSize; i < prevDOFsSize + DOFs.size(); i++) {
						elements[e]->DOFsDomainsCounters()[c * (prevDOFsSize + DOFs.size()) + i] = elements[e]->numberOfLocalDomainsWithDOF(i);
					}
				}
			}

			if (elements[e]->clusters().size() > 1) {
				for (size_t c = 0; c < elements[e]->clusters().size(); c++) {
					if (elements[e]->clusters()[c] == environment->MPIrank) {
						continue;
					}

					sBuffer[t][n2i(elements[e]->clusters()[c])].push_back(elements[e]->vtkCode());
					for (size_t n = 0; n < elements[e]->coarseNodes(); n++) {
						sBuffer[t][n2i(elements[e]->clusters()[c])].push_back(l2g[elements[e]->node(n)]);
					}

					sBuffer[t][n2i(elements[e]->clusters()[c])].push_back(elements[e]->domains().size());
					for (size_t i = prevDOFsSize; i < prevDOFsSize + DOFs.size(); i++) {
						sBuffer[t][n2i(elements[e]->clusters()[c])].push_back(elements[e]->DOFsDomainsCounters()[cluster * (prevDOFsSize + DOFs.size()) + i]);
					}
				}
			}

		}
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < neighbours.size(); n++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, neighbours)) {
		ESINFO(ERROR) << "problem while synchronization of DOFs counters.";
	}

	// Create list of neighbours elements
	std::vector<std::vector<Element*> > nElements(neighbours.size());
	#pragma omp parallel for
	for (size_t n = 0; n < neighbours.size(); n++) {
		size_t p = 0;
		while (p + 1 < rBuffer[n].size()) {
			switch (rBuffer[n][p++]) {
			case NodeVTKCode:
				nElements[n].push_back(new Node(rBuffer[n][p]));
				break;
			case DOFVTKCode:
				nElements[n].push_back(new DOF(rBuffer[n][p]));
				break;
			case Line2VTKCode:
			case Line3VTKCode:
				nElements[n].push_back(new Line2(&rBuffer[n][p]));
				break;
			case Square4VTKCode:
			case Square8VTKCode:
				nElements[n].push_back(new Square4(&rBuffer[n][p]));
				break;
			case Triangle3VTKCode:
			case Triangle6VTKCode:
				nElements[n].push_back(new Triangle3(&rBuffer[n][p]));
				break;
			case UnknownPointVTKCode:
			case UnknownLineVTKCode:
			case UnknownPlaneVTKCode:
				ESINFO(GLOBAL_ERROR) << "Unknown elements are not allowed to send";
				break;
			default:
				// Volume elements are never exchanged
				ESINFO(GLOBAL_ERROR) << "Unknown neighbour element";
			}
			p += nElements[n].back()->nodes();
			for (size_t i = 0; i < nElements[n].back()->nodes(); i++) {
				nElements[n].back()->node(i) = std::lower_bound(g2l.begin(), g2l.end(), nElements[n].back()->node(i), [] (const G2L &mapping, esglobal index) {
					return mapping.global < index;
				})->local;
			}

			nElements[n].back()->numberOfGlobalDomains(rBuffer[n][p++]);
			nElements[n].back()->DOFsDomainsCounters() = std::vector<eslocal>(&rBuffer[n][p], &rBuffer[n][p] + DOFs.size());
			p += DOFs.size();
		}
	}

	// TODO: parallelization
	for (size_t n = 0; n < neighbours.size(); n++) {
		for (size_t e = 0, offset = 0; e < nElements[n].size(); e++) {
			auto it = std::lower_bound(elements.begin(), elements.end(), nElements[n][e], [&] (Element *el1, Element *el2) { return *el1 < *el2; });
			if (it != elements.end() && **it == *(nElements[n][e])) {

				size_t cluster = std::lower_bound((*it)->clusters().begin(), (*it)->clusters().end(), neighbours[n]) - (*it)->clusters().begin();

				(*it)->numberOfGlobalDomains((*it)->numberOfGlobalDomains() + nElements[n][e]->numberOfGlobalDomains());
				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					(*it)->DOFsDomainsCounters()[cluster * (prevDOFsSize + DOFs.size()) + prevDOFsSize + dof] = nElements[n][e]->DOFsDomainsCounters()[dof];
				}
				if ((*it)->clusterOffsets().size() == 0) {
					(*it)->clusterOffsets().resize((*it)->clusters().size());
				}
				(*it)->clusterOffsets()[cluster] = offset++;
			}
		}
	}

	#pragma omp parallel for
	for (size_t n = 0; n < neighbours.size(); n++) {
		for (size_t e = 0; e < nElements[n].size(); e++) {
			delete nElements[n][e];
		}
	}

	// Remove elements that are not in both clusters
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (elements[e]->clusters().size() > 1) {
				std::vector<eslocal> &counters = elements[e]->DOFsDomainsCounters();
				for (size_t c = 0; c < elements[e]->clusters().size(); c++) {
					if (std::all_of(counters.begin() + c * (prevDOFsSize + DOFs.size()), counters.begin() + (c + 1) * (prevDOFsSize + DOFs.size()), [] (eslocal &c) { return c == -1; })) {
						counters.erase(counters.begin() + c * (prevDOFsSize + DOFs.size()), counters.begin() + (c + 1) * (prevDOFsSize + DOFs.size()));
						elements[e]->clusters().erase(elements[e]->clusters().begin() + c--);
					}
				}
			}

		}
	}
}


void Mesh::computeNodesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_nodes, DOFs, _neighbours, _coordinates->_globalIndex, _coordinates->_globalMapping);
}

void Mesh::computeEdgesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_edges, DOFs, _neighbours, _coordinates->_globalIndex, _coordinates->_globalMapping);
}

void Mesh::computeFacesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_faces, DOFs, _neighbours, _coordinates->_globalIndex, _coordinates->_globalMapping);
}

void Mesh::clearNodesDOFsCounters()
{
	for (size_t e = 0; e < _nodes.size(); e++) {
		_nodes[e]->DOFsDomainsCounters().clear();
	}
}

void APIMesh::computeDOFsDOFsCounters()
{
	computeDOFsCounters(_DOFs, { Property::UNKNOWN }, _neighbours, _l2g, *_g2l);
}

void Mesh::mapCoordinatesToDomains()
{
	_coordinates->_clusterIndex.clear();
	_coordinates->_clusterIndex.resize(parts());

	for (size_t p = 0; p < parts(); p++) {
		std::vector<eslocal> l2g;
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			l2g.insert(l2g.end(), _elements[e]->indices(), _elements[e]->indices() + _elements[e]->nodes());
		}

		std::sort(l2g.begin(), l2g.end());
		Esutils::removeDuplicity(l2g);

		_coordinates->_clusterIndex[p] = l2g;
	}
}

struct __Point__ {

	static constexpr size_t digits = 1000;

	__Point__(): x(0), y(0), z(0), id(-1) {};
	__Point__(const Point &p, esglobal id): x(p.x), y(p.y), z(p.z), id(id) {};

	bool operator<(const __Point__ &p) const
	{
		if (std::trunc(z * digits) == std::trunc(p.z * digits)) {
			if (std::trunc(y * digits) == std::trunc(p.y * digits)) {
				if (std::trunc(x * digits) == std::trunc(p.x * digits)) {
					return false;
				}
				return x < p.x;
			}
			return y < p.y;
		}
		return z < p.z;
	}

	bool operator==(const __Point__ &p) const
	{
		return !(*this < p) && !(p < *this);
	}

	double x, y, z;
	esglobal id;
};

void Mesh::synchronizeGlobalIndices()
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_neighbours.begin(), _neighbours.end(), neighbour) - _neighbours.begin();
	};

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _nodes.size());

	// clusters x threads x nodes
	std::vector<std::vector<std::vector<__Point__> > > sBuffer(threads, std::vector<std::vector<__Point__> >(_neighbours.size()));
	std::vector<std::vector<__Point__> > rBuffer(_neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {

			if (_nodes[n]->clusters().size() > 1) {
				if (_nodes[n]->clusters().front() == environment->MPIrank) {
					for (size_t c = 1; c < _nodes[n]->clusters().size(); c++) {
						sBuffer[t][n2i(_nodes[n]->clusters()[c])].push_back(__Point__(_coordinates->_points[n], _coordinates->globalIndex(n)));
					}
				} else {
					sBuffer[t][n2i(_nodes[n]->clusters().front())].push_back(__Point__(_coordinates->_points[n], _coordinates->globalIndex(n)));
				}
			}

		}
	}

	#pragma omp parallel for
	for (size_t n = 0; n < _neighbours.size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
		if (_neighbours[n] < environment->MPIrank) {
			rBuffer[n].resize(sBuffer[0][n].size());
			std::sort(sBuffer[0][n].begin(), sBuffer[0][n].end());
		}
	}

	if (!Communication::receiveLowerKnownSize(sBuffer[0], rBuffer, _neighbours)) {
		ESINFO(ERROR) << "problem while synchronization of global indices.";
	}

	for (size_t n = 0; n < _neighbours.size(); n++) {
		if (_neighbours[n] < environment->MPIrank) {
			for (size_t p = 0; p < rBuffer[n].size(); p++) {
				auto it = std::lower_bound(sBuffer[0][n].begin(), sBuffer[0][n].end(), rBuffer[n][p]);
				if (*it == rBuffer[n][p]) {
					_coordinates->_globalIndex[_coordinates->clusterIndex(it->id)] = rBuffer[n][p].id;
				} else {
					ESINFO(ERROR) << "Internal ERROR while synchronization global indices: " << _neighbours[n] << " on " << environment->MPIrank;
				}
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
			_coordinates->_globalMapping[n].global = _coordinates->_globalIndex[n];
			_coordinates->_globalMapping[n].local = n;
		}
	}

	std::sort(_coordinates->_globalMapping.begin(), _coordinates->_globalMapping.end());
}

void Mesh::synchronizeNeighbours()
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<std::vector<std::vector<int> > > sNeighbours(threads, std::vector<std::vector<int> >(environment->MPIsize));
	std::vector<std::vector<int> > rNeighbours;
	_neighbours.clear();

	std::vector<size_t> distribution = Esutils::getDistribution(threads, _nodes.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {

			for (size_t c1 = 0; c1 < _nodes[n]->clusters().size(); c1++) {
				for (size_t c2 = 0; c2 < _nodes[n]->clusters().size(); c2++) {
					if (c1 != c2 && _nodes[n]->clusters()[c1] != environment->MPIrank) {
						sNeighbours[t][_nodes[n]->clusters()[c1]].push_back(_nodes[n]->clusters()[c2]);
					}
				}
			}
		}
		for (int r = 0; r < environment->MPIsize; r++) {
			std::sort(sNeighbours[t][r].begin(), sNeighbours[t][r].end());
			Esutils::removeDuplicity(sNeighbours[t][r]);
		}
	}

	#pragma omp parallel for
	for (int r = 0; r < environment->MPIsize; r++) {
		for (size_t t = 1; t < threads; t++) {
			sNeighbours[0][r].insert(sNeighbours[0][r].end(), sNeighbours[t][r].begin(), sNeighbours[t][r].end());
		}
		std::sort(sNeighbours[0][r].begin(), sNeighbours[0][r].end());
		Esutils::removeDuplicity(sNeighbours[0][r]);
	}
	for (int r = 0; r < environment->MPIsize; r++) {
		if (sNeighbours[0][r].size()) {
			sNeighbours[0][_neighbours.size()].swap(sNeighbours[0][r]);
			_neighbours.push_back(r);
		}
	}
	rNeighbours.resize(_neighbours.size());

	if (!Communication::exchangeUnknownSize(sNeighbours[0], rNeighbours, _neighbours)) {
		ESINFO(ERROR) << "problem while synchronization of neighbours clusters.";
	}

	size_t size = _neighbours.size();
	for (size_t n = 0; n < size; n++) {
		_neighbours.insert(_neighbours.end(), rNeighbours[n].begin(), rNeighbours[n].end());
	}
	std::sort(_neighbours.begin(), _neighbours.end());
	Esutils::removeDuplicity(_neighbours);

	std::vector<std::vector<std::vector<esglobal> > > sIndices(threads, std::vector<std::vector<esglobal> >(_neighbours.size()));
	std::vector<std::vector<esglobal> > rIndices(_neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
			for (size_t i = 0; i < _neighbours.size(); i++) {
				sIndices[t][i].push_back(_coordinates->globalIndex(n));
			}
		}
	}

	for (size_t n = 0; n < _neighbours.size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sIndices[0][n].insert(sIndices[0][n].end(), sIndices[t][n].begin(), sIndices[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sIndices[0], rIndices, _neighbours)) {
		ESINFO(ERROR) << "problem while synchronization of indices in synchronization neighbours.";
	}

	for (size_t n = 0; n < _neighbours.size(); n++) {
		#pragma omp parallel for
		for (size_t i = 0; i < rIndices[n].size(); i++) {
			eslocal index = _coordinates->clusterIndex(rIndices[n][i]);
			if (index >= 0) {
				Element *node = _nodes[index];
				node->clusters().push_back(_neighbours[n]);
				std::sort(node->clusters().begin(), node->clusters().end());
				Esutils::removeDuplicity(node->clusters());
			}
		}
	}
}

void Mesh::synchronizeRegionOrder()
{
	std::vector<Region*> regions(_regions.begin(), _regions.begin() + 2);
	std::vector<char> name;
	for (size_t r = 2; r < _regions.size(); r++) {
		name = std::vector<char>(_regions[r]->name.begin(), _regions[r]->name.end());
		Communication::broadcastUnknownSize(name);
		regions.push_back(region(std::string(name.begin(), name.end())));
	}
	_regions.swap(regions);
}

void Mesh::checkNeighbours()
{
	ESINFO(ALWAYS) << Info::TextColor::BLUE << "Checking whether neighbours are correct";
	int nSize = _neighbours.size();
	std::vector<int> neighbours;
	std::vector<int> counters(environment->MPIsize);
	std::vector<int> displacements;

	MPI_Gather(&nSize, 1, MPI_INT, counters.data(), 1, MPI_INT, 0, environment->MPICommunicator);
	for (size_t i = 0; i < counters.size(); i++) {
		displacements.push_back(displacements.size() ? displacements.back() + counters[i - 1] : 0);
		neighbours.resize(neighbours.size() + counters[i]);
	}
	MPI_Gatherv(_neighbours.data(), _neighbours.size(), MPI_INT, neighbours.data(), counters.data(), displacements.data(), MPI_INT, 0, environment->MPICommunicator);

	for (int r = 0; r < environment->MPIsize; r++) {
		for (int n = displacements[r]; n < displacements[r] + counters[r]; n++) {
			if (!std::binary_search(neighbours.begin() + displacements[neighbours[n]], neighbours.begin() + displacements[neighbours[n]] + counters[neighbours[n]], r)) {
				ESINFO(GLOBAL_ERROR) << "ESPRESO INTERNAL TEST FAILED: neighbours are not correctly set.";
			}
		}
	}

	std::vector<esglobal> sClusters;
	std::vector<esglobal> rClusters;
	displacements.clear();

	for (size_t n = 0; n < _nodes.size(); n++) {
		sClusters.push_back(_coordinates->globalIndex(n));
		sClusters.push_back(_nodes[n]->clusters().size());
		sClusters.insert(sClusters.end(), _nodes[n]->clusters().begin(), _nodes[n]->clusters().end());
	}
	nSize = sClusters.size();

	MPI_Gather(&nSize, 1, MPI_INT, counters.data(), 1, MPI_INT, 0, environment->MPICommunicator);
	for (size_t i = 0; i < counters.size(); i++) {
		displacements.push_back(displacements.size() ? displacements.back() + counters[i - 1] : 0);
		rClusters.resize(rClusters.size() + counters[i]);
		counters[i] *= sizeof(esglobal);
	}
	MPI_Gatherv(sClusters.data(), sizeof(esglobal) * sClusters.size(), MPI_BYTE, rClusters.data(), counters.data(), displacements.data(), MPI_BYTE, 0, environment->MPICommunicator);

	std::vector<std::vector<int> > nodes;
	for (size_t n = 0; n < rClusters.size(); ) {
		esglobal index = rClusters[n++];
		esglobal cSize = rClusters[n++];
		std::vector<int> clusters;
		for (esglobal c = 0; c < cSize; c++) {
			clusters.push_back(rClusters[n++]);
		}

		if (index + 1 >= (esglobal)nodes.size()) {
			nodes.resize(index + 1);
		}
		if (nodes[index].size()) {
			if (clusters != nodes[index]) {
				ESINFO(GLOBAL_ERROR) << "ESPRESO INTERNAL TEST FAILED: neighbours for node '" << index << "' are not correct.";
			}
		} else {
			nodes[index].swap(clusters);
		}
	}
}

void Mesh::storeNodeData(const std::string &name, std::function<void (std::ofstream &os, const Element* e)> store)
{
	ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing node data: '" << name << "'";
	std::ofstream os(Logging::prepareFile(name));
	for (size_t n = 0; n < _nodes.size(); n++) {
		os << _coordinates->globalIndex(n) << " :: ";
		store(os, _nodes[n]);
	}
	os.close();
}

void Mesh::storeRegions()
{
	ESINFO(ALWAYS) << Info::TextColor::BLUE << "Storing regions";
	for (size_t r = 0; r < _regions.size(); r++) {
		std::ofstream os(Logging::prepareFile(std::to_string(r) + "_" + _regions[r]->name));
		for (size_t e = 0; e < _regions[r]->elements().size(); e++) {
			os << *_regions[r]->elements()[e] << "\n";
		}
		os.close();
	}
}

void Mesh::checkRegions(const std::vector<Element*> &elements)
{
	storeNodeData("clusters", [] (std::ofstream &os, const Element* e) { os << e->clusters(); });
	storeNodeData("regions", [] (std::ofstream &os, const Element* e) {
		std::for_each(e->regions().begin(), e->regions().end(), [&] (const Region *r) { os << r->name << " "; });
		os << "\n";
	});
	storeRegions();
	ESINFO(ALWAYS) << Info::TextColor::BLUE << "Checking whether regions are correct";

	int rSize = _regions.size();
	std::vector<int> sizes(environment->MPIsize);

	MPI_Gather(&rSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, environment->MPICommunicator);
	if (!std::all_of(sizes.begin(), sizes.end(), [&] (int size) { return size == rSize; })) {
		ESINFO(GLOBAL_ERROR) << "ESPRESO INTERNAL ERROR: The same regions have to be set for all processes.";
	}

	std::vector<char> regions;
	for (size_t r = 0; r < _regions.size(); r++) {
		rSize = _regions[r]->name.size();
		MPI_Gather(&rSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, environment->MPICommunicator);
		if (!std::all_of(sizes.begin(), sizes.end(), [&] (int size) { return size == rSize; })) {
			ESINFO(GLOBAL_ERROR) << "ESPRESO INTERNAL ERROR: regions have not the same names on all processes.";
		}
		regions.resize(rSize * environment->MPIsize);
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Gather(const_cast<char*>(_regions[r]->name.c_str()), rSize, MPI_BYTE, regions.data(), rSize, MPI_BYTE, 0, environment->MPICommunicator);
		if (!environment->MPIrank) {
			for (int i = 0; i < environment->MPIsize; i++) {
				std::string str(regions.begin() + i * rSize, regions.begin() + (i + 1) * rSize);
				if (str != _regions[r]->name) {
					ESINFO(GLOBAL_ERROR) << "ESPRESO INTERNAL ERROR: regions have not the same names on all processes.";
				}
			}
		}
	}

	std::vector<esglobal> sRegions;
	std::vector<esglobal> rRegions;
	std::vector<int> displacements;

	auto r2i = [ & ] (const Region *region) -> int {
		for (size_t r = 0; r < _regions.size(); r++) {
			if (_regions[r] == region) {
				return r;
			}
		}
		return -1;
	};

	for (size_t n = 0; n < _nodes.size(); n++) {
		sRegions.push_back(_coordinates->globalIndex(n));
		sRegions.push_back(_nodes[n]->regions().size());
		for (size_t r = 0; r < _nodes[n]->regions().size(); r++) {
			sRegions.push_back(r2i(_nodes[n]->regions()[r]));
		}
	}

	rSize = sRegions.size();

	MPI_Gather(&rSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, environment->MPICommunicator);
	for (size_t i = 0; i < sizes.size(); i++) {
		displacements.push_back(displacements.size() ? displacements.back() + sizes[i - 1] : 0);
		rRegions.resize(rRegions.size() + sizes[i]);
		sizes[i] *= sizeof(esglobal);
	}
	MPI_Gatherv(sRegions.data(), sizeof(esglobal) * sRegions.size(), MPI_BYTE, rRegions.data(), sizes.data(), displacements.data(), MPI_BYTE, 0, environment->MPICommunicator);

	std::vector<std::vector<int> > nodes;
	for (size_t n = 0; n < rRegions.size(); ) {
		esglobal index = rRegions[n++];
		esglobal size = rRegions[n++];
		std::vector<int> regions;
		for (esglobal c = 0; c < size; c++) {
			regions.push_back(rRegions[n++]);
		}
		std::sort(regions.begin(), regions.end());

		if (index + 1 >= (esglobal)nodes.size()) {
			nodes.resize(index + 1);
		}
		if (nodes[index].size()) {
			if (regions != nodes[index]) {
				ESINFO(GLOBAL_ERROR) << "ESPRESO INTERNAL TEST FAILED: regions for node '" << index << "' are not correct.";
			}
		} else {
			nodes[index].swap(regions);
		}
	}

}

}



