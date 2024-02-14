
#include "block.h"
#include "blocksettings.h"

#include "input/parsers/meshgenerator/meshgenerator.h"
#include "input/parsers/meshgenerator/selection/blockborder.h"
#include "input/meshbuilder.h"

#include "input/parsers/meshgenerator/elements/2D/square4.h"
#include "input/parsers/meshgenerator/elements/2D/square8.h"
#include "input/parsers/meshgenerator/elements/2D/triangle3.h"
#include "input/parsers/meshgenerator/elements/2D/triangle6.h"

#include "input/parsers/meshgenerator/elements/3D/tetrahedron4.h"
#include "input/parsers/meshgenerator/elements/3D/tetrahedron10.h"
#include "input/parsers/meshgenerator/elements/3D/pyramid5.h"
#include "input/parsers/meshgenerator/elements/3D/pyramid13.h"
#include "input/parsers/meshgenerator/elements/3D/prisma6.h"
#include "input/parsers/meshgenerator/elements/3D/prisma15.h"
#include "input/parsers/meshgenerator/elements/3D/hexahedron8.h"
#include "input/parsers/meshgenerator/elements/3D/hexahedron20.h"

#include "basis/containers/point.h"
#include "basis/evaluator/evaluator.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "config/ecf/input/block.h"
#include "config/ecf/input/generatorelements.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

BlockGenerator::BlockGenerator(const BlockGeneratorConfiguration &configuration, const BlockSettings &block)
: _block(block)
{
	switch (configuration.element_type) {
	case GENERATOR_ELEMENT_TYPE::SQUARE4:
		_element = new Square4Generator();
		info::mesh->dimension = 2;
		break;
	case GENERATOR_ELEMENT_TYPE::SQUARE8:
		_element = new Square8Generator();
		info::mesh->dimension = 2;
		break;
	case GENERATOR_ELEMENT_TYPE::TRIANGLE3:
		_element = new Triangle3Generator();
		info::mesh->dimension = 2;
		break;
	case GENERATOR_ELEMENT_TYPE::TRIANGLE6:
		_element = new Triangle6Generator();
		info::mesh->dimension = 2;
		break;
	case GENERATOR_ELEMENT_TYPE::TETRA4:
		_element = new Tetrahedron4Generator();
		info::mesh->dimension = 3;
		break;
	case GENERATOR_ELEMENT_TYPE::PYRAMID5:
		_element = new Pyramid5Generator();
		info::mesh->dimension = 3;
		break;
	case GENERATOR_ELEMENT_TYPE::PRISMA6:
		_element = new Prisma6Generator();
		info::mesh->dimension = 3;
		break;
	case GENERATOR_ELEMENT_TYPE::HEXA8:
		_element = new Hexahedron8Generator();
		info::mesh->dimension = 3;
		break;
	case GENERATOR_ELEMENT_TYPE::TETRA10:
		_element = new Tetrahedron10Generator();
		info::mesh->dimension = 3;
		break;
	case GENERATOR_ELEMENT_TYPE::PYRAMID13:
		_element = new Pyramid13Generator();
		info::mesh->dimension = 3;
		break;
	case GENERATOR_ELEMENT_TYPE::PRISMA15:
		_element = new Prisma15Generator();
		info::mesh->dimension = 3;
		break;
	case GENERATOR_ELEMENT_TYPE::HEXA20:
		_element = new Hexahedron20Generator();
		info::mesh->dimension = 3;
		break;
	default:
		eslog::error("Invalid element type.\n");
	}
}

BlockGenerator::~BlockGenerator()
{
	delete _element;
}

void BlockGenerator::coordinates(MeshBuilder &mesh)
{
	Triple<size_t> nodes = _block.domains * _block.elements * (Triple<size_t>(_element->subnodes) - 1) + 1;
	Triple<double> step = (_block.end - _block.start) / Triple<double>(nodes.x, nodes.y, nodes.z).steps();

	mesh.coordinates.reserve(nodes.mul());
	double &xx = _block.projection_x->getX(), &xy = _block.projection_x->getY(), &xz = _block.projection_x->getZ();
	double &yx = _block.projection_y->getX(), &yy = _block.projection_y->getY(), &yz = _block.projection_y->getZ();
	double &zx = _block.projection_z->getX(), &zy = _block.projection_z->getY(), &zz = _block.projection_z->getZ();
	for (size_t z = 0; z < nodes.z; z++) {
		for (size_t y = 0; y < nodes.y; y++) {
			for (size_t x = 0; x < nodes.x; x++) {
				xx = yx = zx = MeshGenerator::precision * (_block.start.x + x * step.x);
				xy = yy = zy = MeshGenerator::precision * (_block.start.y + y * step.y);
				xz = yz = zz = MeshGenerator::precision * (_block.start.z + z * step.z);
				mesh.coordinates.push_back(Point(_block.projection_x->evaluate(), _block.projection_y->evaluate(), _block.projection_z->evaluate()));
			}
		}
	}
}

void BlockGenerator::forEachElement(
		const Triple<size_t> &begin,
		const Triple<size_t> &end,
		std::function<void(std::vector<esint> &indices)> operation,
		std::function<void(Triple<size_t> &offset)> restriction)
{
	Triple<size_t> enodes = Triple<size_t>(_element->subnodes) - 1;
	Triple<size_t> size = (_block.domains * _block.elements * enodes + 1).toSize();

	Triple<size_t> _begin = begin;
	Triple<size_t> _end = end;

	auto correct = [] (size_t &s, size_t &e) {
		if (s == e) { if (s == 0) { e++; } else { s--; } }
	};

	correct(_begin.x, _end.x);
	correct(_begin.y, _end.y);
	correct(_begin.z, _end.z);

	std::vector<esint> indices((enodes + 1).mul());

	Triple<size_t> e;
	for (e.z = _begin.z; e.z < _end.z; e.z++) {
		for (e.y = _begin.y; e.y < _end.y; e.y++) {
			for (e.x = _begin.x; e.x < _end.x; e.x++) {

				size_t index = 0;
				Triple<size_t> eoffset, noffset = e * enodes;
				for (eoffset.z = 0; eoffset.z <= enodes.z; eoffset.z++) {
					for (eoffset.y = 0; eoffset.y <= enodes.y; eoffset.y++) {
						for (eoffset.x = 0; eoffset.x <= enodes.x; eoffset.x++) {
							Triple<size_t> coffset = eoffset + noffset;
							restriction(coffset);
							indices[index++] = (coffset * size).sum();
						}
					}
				}
				operation(indices);
			}
		}
	}
}

void BlockGenerator::elements(MeshBuilder &mesh)
{
	mesh.esize.resize(_element->subelements * (_block.domains * _block.elements).mul(), _element->enodes);
	mesh.etype.resize(_element->subelements * (_block.domains * _block.elements).mul(), (esint)_element->code);
	mesh.enodes.reserve(mesh.etype.size() * _element->enodes);

	auto add = [&] (std::vector<esint> &indices) {
		_element->pushElements(mesh.enodes, indices);
	};

	auto restriction = [] (Triple<size_t> &offset) {};

	Triple<size_t> offset;
	for (offset.z = 0; offset.z < _block.domains.z; offset.z++) {
		for (offset.y = 0; offset.y < _block.domains.y; offset.y++) {
			for (offset.x = 0; offset.x < _block.domains.x; offset.x++) {
				forEachElement(offset * _block.elements, (offset + 1) * _block.elements, add, restriction);
			}
		}
	}
}

void BlockGenerator::neighbors(const std::vector<int> &surroundings, MeshBuilder &mesh)
{
	std::vector<int> sorted(surroundings.size());
	std::iota(sorted.begin(), sorted.end(), 0);
	std::sort(sorted.begin(), sorted.end(), [&] (int i, int j) { return surroundings[i] < surroundings[j]; });

	Triple<size_t> count = _block.domains * _block.elements * (Triple<size_t>(_element->subnodes) - 1);
	Triple<size_t> size = (count + 1).toSize();

	std::vector<std::pair<esint, esint> > neighs;

	for (int i = 0; i < 27; i++) {
		if (surroundings[sorted[i]] < 0) {
			continue;
		}
		size_t  sx = sorted[i] % 3 == 2 ? count.x : 0,
				ex = sorted[i] % 3 != 0 ? count.x : 0,
				sy = sorted[i] % 9 >= 6 ? count.y : 0,
				ey = sorted[i] % 9 >= 3 ? count.y : 0,
				sz = sorted[i] / 9 == 2 ? count.z : 0,
				ez = sorted[i] / 9 >= 1 ? count.z : 0;


		Triple<size_t> offset;
		for (offset.x = sx; offset.x <= ex; offset.x++) {
			for (offset.y = sy; offset.y <= ey; offset.y++) {
				for (offset.z = sz; offset.z <= ez; offset.z++) {
					neighs.push_back(std::make_pair((offset * size).sum(), surroundings[sorted[i]]));
				}
			}
		}
	}

	std::sort(neighs.begin(), neighs.end());

	mesh._nrankdist.reserve(mesh.nIDs.size() + 1);
	mesh._nrankdist.push_back(0);
	mesh._nranks.reserve(neighs.size());
	for (size_t i = 0; i < neighs.size(); i++) {
		if (i && neighs[i - 1].first != neighs[i].first) {
			mesh._nrankdist.push_back(mesh._nranks.size());
		}
		mesh._nranks.push_back(neighs[i].second);
	}
	mesh._nrankdist.push_back(mesh._nranks.size());
}

bool BlockGenerator::region(
		BlockBorder &intersection,
		Triple<size_t> &ebegin, Triple<size_t> &eend,
		Triple<size_t> &nbegin, Triple<size_t> &nend)
{
	eend = _block.domains * _block.elements;
	nend = _block.domains * _block.elements * (Triple<size_t>(_element->subnodes) - 1);
	Triple<size_t> lastNode = nend;

	// project intersection to interval <0, 1>
	Triple<double> istart = (Triple<double>)(intersection.start - _block.start) / (_block.end - _block.start);
	Triple<double> iend   = (Triple<double>)(intersection.end - _block.start) / (_block.end - _block.start);

	// get start, end node
	Triple<double> nStart = istart * (nend - nbegin);
	Triple<double> nEnd   = iend   * (nend - nbegin);
	// shrink interval to the closest nodes
	nbegin = nStart.ceil();
	nend = nEnd.floor();

	// If interval is open, checks whether min, max node should be modev
	Triple<double> cStart = nStart.ceil();
	Triple<double> fEnd   = nEnd.floor();
	nStart = (nStart / MeshGenerator::precision).ceil();
	nEnd   = (nEnd   / MeshGenerator::precision).floor();

	if (intersection.excludeStart.x && cStart.x / MeshGenerator::precision == nStart.x) { nbegin.x++; }
	if (intersection.excludeStart.y && cStart.y / MeshGenerator::precision == nStart.y) { nbegin.y++; }
	if (intersection.excludeStart.z && cStart.z / MeshGenerator::precision == nStart.z) { nbegin.z++; }
	if (intersection.excludeEnd.x && fEnd.x / MeshGenerator::precision == nEnd.x) { nend.x--; }
	if (intersection.excludeEnd.y && fEnd.y / MeshGenerator::precision == nEnd.y) { nend.y--; }
	if (intersection.excludeEnd.z && fEnd.z / MeshGenerator::precision == nEnd.z) { nend.z--; }

	// No matter what above happen, neighbor nodes have to match! Check it again
	if (intersection.start.x == _block.start.x) { nbegin.x = intersection.excludeStart.x ? 1 : 0; }
	if (intersection.start.y == _block.start.y) { nbegin.y = intersection.excludeStart.y ? 1 : 0; }
	if (intersection.start.z == _block.start.z) { nbegin.z = intersection.excludeStart.z ? 1 : 0; }
	if (intersection.end.x == _block.end.x) { nend.x = intersection.excludeEnd.x ? lastNode.x - 1 : lastNode.x; }
	if (intersection.end.y == _block.end.y) { nend.y = intersection.excludeEnd.y ? lastNode.y - 1 : lastNode.y; }
	if (intersection.end.z == _block.end.z) { nend.z = intersection.excludeEnd.z ? lastNode.z - 1 : lastNode.z; }

	Triple<size_t> enodes = Triple<size_t>(_element->subnodes);
	if (enodes.z == 1) { enodes.z++; }

	// Add only elements that have included more that one node
	ebegin = nbegin / (enodes - 1);
	eend = ((Triple<double>)(nend) / (enodes - 1)).ceil();
	if (_element->subnodes[2] == 1) {
		eend.z++;
	}

	// There is still change to miss all nodes by interval.
	if (nbegin.x > nend.x) { return false; }
	if (nbegin.y > nend.y) { return false; }
	if (nbegin.z > nend.z) { return false; }

	return true;
}

void BlockGenerator::nodesRegion(const BlockBorder &border, std::vector<esint> &nodes)
{
	BlockBorder intersection = border.intersect(_block);
	Triple<size_t> ebegin, eend, nbegin, nend;
	if (!region(intersection, ebegin, eend, nbegin, nend)) {
		return;
	}

	CubeEdge edge;
	auto pushEdge = [&] (std::vector<esint> &indices) {
		_element->pushNodes(nodes, indices, edge);
	};

	CubeFace face;
	auto pushFace = [&] (std::vector<esint> &indices) {
		_element->pushNodes(nodes, indices, face);
	};

	auto restriction = [&] (Triple<size_t> &offset) {
		offset.x = offset.x < nbegin.x ? nbegin.x : nend.x < offset.x ? nend.x : offset.x;
		offset.y = offset.y < nbegin.y ? nbegin.y : nend.y < offset.y ? nend.y : offset.y;
		offset.z = offset.z < nbegin.z ? nbegin.z : nend.z < offset.z ? nend.z : offset.z;
	};

	switch (intersection.dimension()) {

	case 1:
		edge = intersection.getEdge(_block);
		forEachElement(ebegin, eend, pushEdge, restriction);
		utils::sortAndRemoveDuplicates(nodes);
		break;
	case 2:
		face = intersection.getFace(_block);
		forEachElement(ebegin, eend, pushFace, restriction);
		utils::sortAndRemoveDuplicates(nodes);
		break;

	case 0:
		// NOT INTERSECTED IN THIS PROCESS
		break;
	case 3:
		eslog::error("Implement a selection of nodes inside 3D interval.\n");
		break;
	default:
		break;
	}
}

void BlockGenerator::edgesRegion(const BlockBorder &border, MeshBuilder &mesh, std::vector<esint> &elements)
{
	BlockBorder intersection = border.intersect(_block);
	Triple<size_t> ebegin, eend, nbegin, nend;
	if (!region(intersection, ebegin, eend, nbegin, nend) || intersection.dimension() < 1) {
		return;
	}

	CubeEdge edge = intersection.getEdge(_block);
	auto pushEdge = [&] (std::vector<esint> &indices) {
		_element->pushEdge(mesh.enodes, mesh.esize, mesh.etype, indices, edge);
	};

	auto restriction = [&] (Triple<size_t> &offset) {};

	size_t begin = mesh.etype.size();

	switch (intersection.dimension()) {

	case 1:
		forEachElement(ebegin, eend, pushEdge, restriction);
		break;

	case 0:
		// NOT INTERSECTED IN THIS PROCESS
		break;
	case 2:
		eslog::error("Implement a selection of edges inside 2D interval.\n");
		break;
	case 3:
		eslog::error("Implement a selection of edges inside 3D interval.\n");
		break;
	}

	elements.resize(mesh.etype.size() - begin);
	std::iota(elements.begin(), elements.end(), begin);
}

void BlockGenerator::facesRegion(const BlockBorder &border, MeshBuilder &mesh, std::vector<esint> &elements)
{
	BlockBorder intersection = border.intersect(_block);
	Triple<size_t> ebegin, eend, nbegin, nend;
	if (!region(intersection, ebegin, eend, nbegin, nend)) {
		return;
	}

	CubeFace face = intersection.getFace(_block);
	auto pushFace = [&] (std::vector<esint> &indices) {
		_element->pushFace(mesh.enodes, mesh.esize, mesh.etype, indices, face);
	};

	auto restriction = [&] (Triple<size_t> &offset) {};

	size_t begin = mesh.etype.size();

	switch (intersection.dimension()) {

	case 2:
		forEachElement(ebegin, eend, pushFace, restriction);
		break;

	case 0:
	case 1:
		// NOT INTERSECTED IN THIS PROCESS
		break;
	case 3:
		eslog::error("Implement a selection of faces inside 3D interval.\n");
		break;
	}

	elements.resize(mesh.etype.size() - begin);
	std::iota(elements.begin(), elements.end(), begin);
}

void BlockGenerator::elementsRegion(const BlockBorder &border, std::vector<esint> &elements)
{
	BlockBorder intersection = border.intersect(_block);
	Triple<size_t> ebegin, eend, nbegin, nend;
	if (!region(intersection, ebegin, eend, nbegin, nend)) {
		return;
	}

	auto isIn = [&] (const Triple<size_t> &index) {
		return
				ebegin.x <= index.x && index.x < eend.x &&
				ebegin.y <= index.y && index.y < eend.y &&
				ebegin.z <= index.z && index.z < eend.z;
	};

	Triple<size_t> domain, element; size_t eindex = 0;
	switch (intersection.dimension()) {

	case 0:
	case 1:
		// NOT INTERSECTED IN THIS PROCESS
		break;
	case 2:
	case 3:
		for (domain.z = 0; domain.z < _block.domains.z; domain.z++) {
			for (domain.y = 0; domain.y < _block.domains.y; domain.y++) {
				for (domain.x = 0; domain.x < _block.domains.x; domain.x++) {

					for (element.z = 0; element.z < _block.elements.z; element.z++) {
						for (element.y = 0; element.y < _block.elements.y; element.y++) {
							for (element.x = 0; element.x < _block.elements.x; element.x++, eindex += _element->subelements) {

								if (isIn(domain * _block.elements + element)) {
									for (size_t e = 0; e < _element->subelements; e++) {
										elements.push_back(eindex + e);
									}
								}

							}
						}
					}

				}
			}
		}
		break;
	default:
		break;
	}
}

void BlockGenerator::pattern(const Triple<size_t> &offset, const Triple<size_t> &size, std::vector<esint> &elements, Pattern pattern, size_t psize)
{
	Triple<size_t> divisor(size * _block.elements * _block.domains / psize);
	Triple<size_t> modulo;

	divisor.z = divisor.z ? divisor.z : 1;

	Triple<size_t> domain, element, begin, index;
	begin = offset * _block.elements * _block.domains;

	size_t color = 0;
	switch (pattern) {
	case Pattern::CHESSBOARD_BLACK:
		color = 1;
		break;
	case Pattern::CHESSBOARD_WHITE:
		color = 0;
		break;
	}
	size_t eindex = 0;
	for (domain.z = 0; domain.z < _block.domains.z; domain.z++) {
		for (domain.y = 0; domain.y < _block.domains.y; domain.y++) {
			for (domain.x = 0; domain.x < _block.domains.x; domain.x++) {

				for (element.z = 0; element.z < _block.elements.z; element.z++) {
					for (element.y = 0; element.y < _block.elements.y; element.y++) {
						for (element.x = 0; element.x < _block.elements.x; element.x++, eindex += _element->subelements) {

							index = begin + domain * _block.elements + element;
							if ((index / divisor).sum() % 2 == color) {
								for (size_t e = 0; e < _element->subelements; e++) {
									elements.push_back(eindex + e);
								}
							}

						}
					}
				}

			}
		}
	}
}


