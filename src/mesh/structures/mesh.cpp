
#include "mesh.h"
#include "mkl.h"
#include "../../config/configuration.h"

namespace espreso {

Mesh::Mesh():_elements(0)
{
	_partPtrs.resize(2);
	_partPtrs[0] = 0;
	_partPtrs[1] = 0;

	_regions.push_back(new Region(_elements));
	_regions.back()->name = "ALL_ELEMENTS";
	_regions.push_back(new Region(_nodes));
	_regions.back()->name = "ALL_NODES";
}

void Mesh::partitiate(size_t parts)
{
	if (parts == 1 && this->parts() == 1) {
		_partPtrs = { 0, (eslocal)_elements.size() };
		mapElementsToDomains();
	} else {
		std::vector<eslocal> ePartition = getPartition(0, _elements.size(), parts);

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
	mapFacesToDomains();
	mapEdgesToDomains();
	mapNodesToDomains();
	mapCoordinatesToDomains();
}

void APIMesh::partitiate(size_t parts)
{
	if (parts == 1 && this->parts() == 1) {
		_partPtrs = { 0, (eslocal)_elements.size() };
		mapElementsToDomains();
	} else {
		std::vector<eslocal> ePartition = getPartition(0, _elements.size(), parts);

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

	return computePartition(e, n, end - begin, _coordinates.clusterSize(), ncommon, parts);
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
	std::vector<std::vector<eslocal> > neighbours(_coordinates.localSize(part));
	for (eslocal i = begin; i < end; i++) {
		if (ePartition[i - begin] == subpart) {
			for (size_t j = 0; j < _elements[i]->nodes(); j++) {
				std::vector<eslocal> neigh = _elements[i]->getNeighbours(j);
				for (size_t k = 0; k < neigh.size(); k++) {
					if (_elements[i]->node(j) < neigh[k]) {
						neighbours[_coordinates.localIndex(_elements[i]->node(j), part)].push_back(_coordinates.localIndex(neigh[k], part));
					}
				}
			}
		}
	}

	return _coordinates.clusterIndex(computeCenter(neighbours), part);
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
}

static bool isOuterFace(
		std::vector<std::vector<eslocal> > &nodesElements,
		std::vector<eslocal> &face)
{
	std::vector<eslocal> result(nodesElements[face[0]]);
	std::vector<eslocal>::iterator it = result.end();

	for (size_t i = 1; i < face.size(); i++) {
		std::vector<eslocal> tmp(result.begin(), it);
		it = std::set_intersection(tmp.begin(), tmp.end(),
				nodesElements[face[i]].begin(), nodesElements[face[i]].end(),
				result.begin());
		if (it - result.begin() == 1) {
			return true;
		}
	}

	return false;
}

void Mesh::getSurface(Mesh &surface) const
{
	// vector of faces in all parts
	std::vector<std::vector<std::vector<eslocal> > > faces(parts());
	// number of elements in all parts
	std::vector<size_t> elementsCount(parts(), 0);

	if (parts() < 1) {
		ESINFO(ERROR) << "Internal error: _partPtrs.size().";
	}

	#pragma omp parallel for
	for  (size_t i = 0; i < parts(); i++) {
		// Compute nodes' adjacent elements
		std::vector<std::vector<eslocal> > nodesElements(_coordinates.localSize(i));
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->nodes(); k++) {
				nodesElements[_elements[j]->node(k)].push_back(j);
			}
		}

		// compute number of elements and fill used nodes
		for (eslocal j = _partPtrs[i]; j < _partPtrs[i + 1]; j++) {
			for (size_t k = 0; k < _elements[j]->faces(); k++) {
				std::vector<eslocal> face(_elements[j]->face(k)->indices(), _elements[j]->face(k)->indices() + _elements[j]->face(k)->nodes());
				if (isOuterFace(nodesElements, face)) {
					for (size_t f = 0; f < face.size(); f++) {
						face[f] = _coordinates.clusterIndex(face[f], i);
					}
					faces[i].push_back(face);
					if (face.size() == 3) {
						elementsCount[i] += 1;
					}
					if (face.size() == 4) {
						elementsCount[i] += 2;
					}
				}
			}
		}
	}

	surface._coordinates = _coordinates;

	size_t count = 0;
	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		count += elementsCount[i];
	}

	surface._elements.reserve(count);
	surface._partPtrs.clear();
	surface._partPtrs.reserve(_partPtrs.size());
	eslocal params[6] = {0, 0, 0, 0, 0, 0};

	// create surface mesh
	surface._partPtrs.push_back(0); //(surface._elements.size());
	for (size_t i = 0; i + 1 < _partPtrs.size(); i++) {
		for (size_t j = 0; j < faces[i].size(); j++) {
			std::vector<eslocal> &face = faces[i][j];
			if (face.size() == 3) {
				surface._elements.push_back(new Triangle3(&face[0], params));
			}
			// divide square to triangles
			if (face.size() == 4) {
				size_t min = 0;
				for (size_t p = 1; p < 4; p++) {
					if (_coordinates[face[p]] < _coordinates[face[min]]) {
						min = p;
					}
				}
				if (min % 2 == 0) {
					surface._elements.push_back(new Triangle3(&face[0], params));
					face[1] = face[0];
					surface._elements.push_back(new Triangle3(&face[1], params));
				} else {
					surface._elements.push_back(new Triangle3(&face[1], params));
					face[2] = face[3];
					surface._elements.push_back(new Triangle3(&face[0], params));
				}
			}
		}
		surface._partPtrs.push_back(surface._elements.size());
	}

	surface.mapCoordinatesToDomains();
	surface._neighbours = _neighbours;
	surface._materials = _materials;
	for (size_t i = 0; i < _evaluators.size(); i++) {
		surface._evaluators.push_back(_evaluators[i]->copy());
	}
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

static void _loadProperty(Mesh &mesh, size_t loadStep, std::vector<Evaluator*> &evaluators, const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, bool distributeToNodes)
{
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

	for (auto it = regions.begin(); it != regions.end(); ++it) {
		Region *region = mesh.region(it->first);
		region->settings.resize(loadStep + 1);
		std::vector<std::string> values = Parser::split(it->second, ",");

		for (size_t p = 0; p < properties.size(); p++) {
			std::string value = properties.size() == 1 ? values[0] : getValue(values, parameters[p]);
			if (!value.size()) {
				continue;
			}

			if (!StringCompare::contains(value, "xyzt")) {
				Expression expr(value, {});
				evaluators.push_back(new ConstEvaluator(expr.evaluate({}), properties[p]));
			} else {
				evaluators.push_back(new CoordinatesEvaluator(value, mesh.coordinates(), properties[p]));
			}

			if (distributeToNodes && region->elements().size() && region->elements()[0]->nodes() > 1) {
				ESINFO(OVERVIEW) << "Set " << properties[p] << " to '" << value << "' for LOAD STEP " << loadStep + 1 << " for nodes of region '" << region->name << "'";
				std::vector<Element*> nodes;
				for (size_t i = 0; i < region->elements().size(); i++) {
					for (size_t n = 0; n < region->elements()[i]->nodes(); n++) {
						nodes.push_back(mesh.nodes()[region->elements()[i]->node(n)]);
					}
				}
				std::sort(nodes.begin(), nodes.end());
				Esutils::removeDuplicity(nodes);

				distribute(nodes, properties[p], evaluators.back(), region);
			} else {
				ESINFO(OVERVIEW) << "Set " << properties[p] << " to '" << value << "' for LOAD STEP " << loadStep + 1 << " for region '" << region->name << "'";
				distribute(region->elements(), properties[p], evaluators.back(), region);
			}
			region->settings[loadStep][properties[p]].push_back(evaluators.back());
		}
	}
}

void Mesh::loadProperty(const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, size_t loadStep) {

	_loadProperty(*this, loadStep, _evaluators, regions, parameters, properties, false);
}

void Mesh::loadNodeProperty(const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties, size_t loadStep)
{
	_loadProperty(*this, loadStep, _evaluators, regions, parameters, properties, true);
}

void Mesh::loadProperty(const ConfigurationVectorMap<std::string, std::string> &property, const std::vector<std::string> &parameters, const std::vector<Property> &properties)
{
	for (auto it = property.configurations.begin(); it != property.configurations.end(); ++it) {
		std::stringstream ss(it->first);
		size_t step;
		ss >> step;
		_loadProperty(*this, step - 1, _evaluators, it->second->values, parameters, properties, false);
	}
}

void Mesh::loadNodeProperty(const ConfigurationVectorMap<std::string, std::string> &property, const std::vector<std::string> &parameters, const std::vector<Property> &properties)
{
	for (auto it = property.configurations.begin(); it != property.configurations.end(); ++it) {
		std::stringstream ss(it->first);
		size_t step;
		ss >> step;
		_loadProperty(*this, step - 1, _evaluators, it->second->values, parameters, properties, true);
	}
}

void Mesh::removeDuplicateRegions()
{
	size_t threads = environment->OMP_NUM_THREADS;

	auto remove = [&] (std::vector<Element*> &elements) {
		std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());
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
			merge(*last, *it);
			delete *it;
		}
	}
	elements.resize(++last - elements.begin());
}

template<typename MergeFunction>
static std::vector<Element*> mergeElements(size_t threads, std::vector<size_t> &distribution, std::vector<std::vector<Element*> > &elements, MergeFunction merge)
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

	std::vector<Element*> created = mergeElements(threads, distribution, edges, [] (Element* e1, Element *e2) {
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

	std::vector<Element*> created = mergeElements(threads, distribution, faces, [] (Element* e1, Element *e2) {
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
	_nodes.reserve(_coordinates.clusterSize());
	for (size_t i = 0; i < _coordinates.clusterSize(); i++) {
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

	_edges = mergeElements(threads, distribution, edges, [] (Element* e1, Element *e2) {
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
				}
			}
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
		const std::vector<Property> &DOFs,
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
					if (elements[i]->DOFsIndices()[d * DOFs.size() + dof] == __HAS_DOF__) {

						elements[i]->DOFsIndices()[d * DOFs.size() + dof] = counters[elements[i]->domains()[d]]++;
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

std::vector<size_t> Mesh::assignVariousDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
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


	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _nodes.size());

	// domains x DOFs x (threads + 1)
	std::vector<std::vector<std::vector<size_t> > > threadsOffsets(parts(), std::vector<std::vector<size_t> >(DOFs.size(), std::vector<size_t>(threads + 1)));

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::vector<size_t> > threadOffset(parts(), std::vector<size_t>(DOFs.size(), 0));
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {
			_nodes[i]->DOFsIndices().resize(DOFs.size() * _nodes[i]->domains().size(), __WITHOUT_DOF__);

			for (size_t d = 0; d < _nodes[i]->domains().size(); d++) {
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

	setDOFsIndices(_nodes, _DOFtoElement, parts(), DOFs, offsets, threadsOffsets);

	return sizes;
}


static std::vector<size_t> fillUniformDOFs(
		std::vector<Element*> &elements,
		std::vector<std::vector<Element*> > &DOFtoElement,
		size_t parts,
		const std::vector<Property> &DOFs,
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
			elements[i]->DOFsIndices().resize(DOFs.size() * elements[i]->domains().size(), 1);

			for (size_t d = 0; d < elements[i]->domains().size(); d++) {
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

	setDOFsIndices(elements, DOFtoElement, parts, DOFs, offsets, threadsOffsets);

	return sizes;
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToNodes(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_nodes, _DOFtoElement, parts(), DOFs, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToEdges(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_edges, _DOFtoElement, parts(), DOFs, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToFaces(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_faces, _DOFtoElement, parts(), DOFs, offsets);
}

std::vector<size_t> Mesh::assignUniformDOFsIndicesToElements(const std::vector<size_t> &offsets, const std::vector<Property> &DOFs)
{
	return fillUniformDOFs(_elements, _DOFtoElement, parts(), DOFs, offsets);
}

std::vector<size_t> APIMesh::distributeDOFsToDomains(const std::vector<size_t> &offsets)
{
	return fillUniformDOFs(_DOFs, _DOFtoElement, parts(), { Property::UNKNOWN }, offsets);
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

	// Compute send buffers
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			elements[e]->DOFsDomainsCounters().resize(DOFs.size() * elements[e]->clusters().size(), -1);
			size_t cluster = std::lower_bound(elements[e]->clusters().begin(), elements[e]->clusters().end(), environment->MPIrank) - elements[e]->clusters().begin();
			for (size_t i = 0; i < DOFs.size(); i++) {
				elements[e]->DOFsDomainsCounters()[cluster * DOFs.size() + i] = elements[e]->numberOfLocalDomainsWithDOF(i);
			}
			if (elements[e]->clusters().size() > 1) {
				for (auto c = elements[e]->clusters().begin(); c != elements[e]->clusters().end(); ++c) {
					if (*c == environment->MPIrank) {
						continue;
					}

					sBuffer[t][n2i(*c)].push_back(elements[e]->vtkCode());
					for (size_t n = 0; n < elements[e]->coarseNodes(); n++) {
						sBuffer[t][n2i(*c)].push_back(l2g[elements[e]->node(n)]);
					}

					for (size_t i = 0; i < DOFs.size(); i++) {
						sBuffer[t][n2i(*c)].push_back(elements[e]->DOFsDomainsCounters()[cluster * DOFs.size() + i]);
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

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		MPI_Isend(sBuffer[0][n].data(), sizeof(esglobal) * sBuffer[0][n].size(), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + n);
	}

	int flag;
	size_t counter = 0;
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(esglobal));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	// All data are exchanged

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

			nElements[n].back()->DOFsDomainsCounters() = std::vector<eslocal>(&rBuffer[n][p], &rBuffer[n][p] + DOFs.size());
			p += DOFs.size();
		}
	}

	// TODO: parallelization
	for (size_t n = 0; n < neighbours.size(); n++) {
		for (size_t e = 0, offset = 0; e < nElements[n].size(); e++) {
			auto it = std::lower_bound(elements.begin(), elements.end(), nElements[n][e], [&] (Element *el1, Element *el2) { return *el1 < *el2; });
			if (it != elements.end() && **it == *(nElements[n][e])) {
				(*it)->clusterOffsets().resize((*it)->clusters().size());
				size_t cluster = std::lower_bound((*it)->clusters().begin(), (*it)->clusters().end(), neighbours[n]) - (*it)->clusters().begin();
				for (size_t dof = 0; dof < DOFs.size(); dof++) {
					(*it)->DOFsDomainsCounters()[cluster * DOFs.size() + dof] = nElements[n][e]->DOFsDomainsCounters()[dof];
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
					if (std::all_of(counters.begin() + c * DOFs.size(), counters.begin() + (c + 1) * DOFs.size(), [] (eslocal &c) { return c == -1; })) {
						counters.erase(counters.begin() + c * DOFs.size(), counters.begin() + (c + 1) * DOFs.size());
						elements[e]->clusters().erase(elements[e]->clusters().begin() + c--);
					}
				}
			}

		}
	}
}


void Mesh::computeNodesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_nodes, DOFs, _neighbours, _coordinates._globalIndex, _coordinates._globalMapping);
}

void Mesh::computeEdgesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_edges, DOFs, _neighbours, _coordinates._globalIndex, _coordinates._globalMapping);
}

void Mesh::computeFacesDOFsCounters(const std::vector<Property> &DOFs)
{
	computeDOFsCounters(_faces, DOFs, _neighbours, _coordinates._globalIndex, _coordinates._globalMapping);
}

void APIMesh::computeDOFsDOFsCounters()
{
	computeDOFsCounters(_DOFs, { Property::UNKNOWN }, _neighbours, _l2g, _g2l);
}

void Mesh::mapCoordinatesToDomains()
{
	_coordinates._clusterIndex.clear();
	_coordinates._clusterIndex.resize(parts());

	for (size_t p = 0; p < parts(); p++) {
		std::vector<eslocal> l2g;
		for (eslocal e = _partPtrs[p]; e < _partPtrs[p + 1]; e++) {
			l2g.insert(l2g.end(), _elements[e]->indices(), _elements[e]->indices() + _elements[e]->nodes());
		}

		std::sort(l2g.begin(), l2g.end());
		Esutils::removeDuplicity(l2g);

		_coordinates._clusterIndex[p] = l2g;
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
	std::vector<std::vector<std::vector<__Point__> > > sBuffer(_neighbours.size(), std::vector<std::vector<__Point__> >(threads));
	std::vector<std::vector<__Point__> > rBuffer(_neighbours.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {

			if (_nodes[n]->clusters().size() > 1) {
				if (_nodes[n]->clusters().front() == environment->MPIrank) {
					for (size_t c = 1; c < _nodes[n]->clusters().size(); c++) {
						sBuffer[n2i(_nodes[n]->clusters()[c])][t].push_back(__Point__(_coordinates[n], _coordinates.globalIndex(n)));
					}
				} else {
					sBuffer[n2i(_nodes[n]->clusters().front())][t].push_back(__Point__(_coordinates[n], _coordinates.globalIndex(n)));
				}
			}

		}
	}

	#pragma omp parallel for
	for (size_t n = 0; n < _neighbours.size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sBuffer[n][0].insert(sBuffer[n][0].end(), sBuffer[n][t].begin(), sBuffer[n][t].end());
		}
		if (_neighbours[n] < environment->MPIrank) {
			rBuffer[n].resize(sBuffer[n][0].size());
			std::sort(sBuffer[n][0].begin(), sBuffer[n][0].end());
		}
	}

	std::vector<MPI_Request> req(_neighbours.size());
	for (size_t n = 0; n < _neighbours.size(); n++) {
		if (_neighbours[n] > environment->MPIrank) {
			MPI_Isend(sBuffer[n][0].data(), sizeof(__Point__) * sBuffer[n][0].size(), MPI_BYTE, _neighbours[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
		if (_neighbours[n] < environment->MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(__Point__) * rBuffer[n].size(), MPI_BYTE, _neighbours[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
	}

	MPI_Waitall(_neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	for (size_t n = 0; n < _neighbours.size(); n++) {
		if (_neighbours[n] < environment->MPIrank) {
			for (size_t p = 0; p < rBuffer[n].size(); p++) {
				auto it = std::lower_bound(sBuffer[n][0].begin(), sBuffer[n][0].end(), rBuffer[n][p]);
				if (*it == rBuffer[n][p]) {
					_coordinates._globalIndex[_coordinates.clusterIndex(it->id)] = rBuffer[n][p].id;
				} else {
					ESINFO(ERROR) << "Internal ERROR while synchronization global indices: " << _neighbours[n] << " on " << environment->MPIrank;
				}
			}
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
			_coordinates._globalMapping[n].global = _coordinates._globalIndex[n];
			_coordinates._globalMapping[n].local = n;
		}
	}

	std::sort(_coordinates._globalMapping.begin(), _coordinates._globalMapping.end());
}

}


