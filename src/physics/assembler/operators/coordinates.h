
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"

#include "gausspoints.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct CoordinatesToElementNodes: public Operator {
	CoordinatesToElementNodes(serializededata<esint, esint>::const_iterator procNodes, ParameterData &ncoordinates, int interval)
	: Operator(interval, ncoordinates.isconst[interval], Link(interval).inputs(1).outputs(ncoordinates)),
	  procNodes(procNodes),
	  ncoordinates(ncoordinates, interval, ndim)
	{

	}

	serializededata<esint, esint>::const_iterator procNodes;
	OutputParameterIterator ncoordinates;

	void operator++()
	{
		++procNodes;
	}
};

struct Coordinates2DToElementNodes: CoordinatesToElementNodes {
	GET_NAME(Coordinates2DToElementNodes)
	using CoordinatesToElementNodes::CoordinatesToElementNodes;

	void operator()()
	{
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
			ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
			ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
		}
	}
};

struct Coordinates3DToElementNodes: CoordinatesToElementNodes {
	GET_NAME(Coordinates3DToElementNodes)
	using CoordinatesToElementNodes::CoordinatesToElementNodes;

	void operator()()
	{
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
			ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
			ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
			ncoordinates[2] = info::mesh->nodes->coordinates->datatarray()[*n].z;
		}
	}
};

template <class Kernel>
struct ElementCoordinates: public ElementOperatorBuilder {
	GET_NAME(ElementCoordinates)

	Kernel &kernel;

	ElementCoordinates(Kernel &kernel): kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		kernel.coords.node.addInputs(info::mesh->nodes->coordinates);
		kernel.coords.gp.addInputs(kernel.coords.node);
		return true;
	}

	void apply(int interval)
	{
		auto procNodes = info::mesh->elements->procNodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
		if (info::mesh->dimension == 2) {
			iterate_elements(Coordinates2DToElementNodes(procNodes, kernel.coords.node, interval));
			iterate_elements_gps<Kernel>(FromNodesToGaussPoints<2>(kernel.integration.N, kernel.coords.node, kernel.coords.gp, interval));
		}
		if (info::mesh->dimension == 3) {
			iterate_elements(Coordinates3DToElementNodes(procNodes, kernel.coords.node, interval));
			iterate_elements_gps<Kernel>(FromNodesToGaussPoints<3>(kernel.integration.N, kernel.coords.node, kernel.coords.gp, interval));
		}
	}
};

template <class Kernel>
struct BoundaryCoordinates: public BoundaryOperatorBuilder {
	GET_NAME(BoundaryCoordinates)

	Kernel &kernel;

	BoundaryCoordinates(Kernel &kernel): kernel(kernel)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			kernel.coords.boundary.node.regions[r].addInputs(info::mesh->nodes->coordinates);
			kernel.coords.boundary.gp.regions[r].addInputs(kernel.coords.boundary.node.regions[r]);
			kernel.coords.boundary.node.regions[r].isset = true;
			kernel.coords.boundary.gp.regions[r].isset = true;
		}
		return true;
	}

	void apply(int region, int interval)
	{
		auto procNodes = info::mesh->boundaryRegions[region]->procNodes->cbegin() + info::mesh->boundaryRegions[region]->eintervals[interval].begin;
		if (info::mesh->dimension == 2) {
			iterate_boundary(Coordinates2DToElementNodes(procNodes, kernel.coords.boundary.node.regions[region], interval), region);
			iterate_boundary_gps<Kernel>(FromNodesToGaussPoints<2>(kernel.integration.boundary.N.regions[region], kernel.coords.boundary.node.regions[region], kernel.coords.boundary.gp.regions[region], interval), region);
		}
		if (info::mesh->dimension == 3) {
			iterate_boundary(Coordinates3DToElementNodes(procNodes, kernel.coords.boundary.node.regions[region], interval), region);
			iterate_boundary_gps<Kernel>(FromNodesToGaussPoints<3>(kernel.integration.boundary.N.regions[region], kernel.coords.boundary.node.regions[region], kernel.coords.boundary.gp.regions[region], interval), region);
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_ */
