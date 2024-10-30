
#include "meshbuilder.h"
#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/evaluator/evaluator.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "input/builders/sequentialinput.h"
#include "input/builders/scatteredinput.h"
#include "input/builders/generatedinput.h"

#include <algorithm>
#include <cstdlib>

//#include "wrappers/mpi/communication.h"
//#include "basis/utilities/print.h"
//#include <fstream>

using namespace espreso;

void MeshBuilder::build()
{
	profiler::syncstart("build_mesh");
	profiler::syncparam("nodes", nIDs.size());
	profiler::syncparam("elements", eIDs.size());
	profiler::syncparam("enodes", enodes.size());
	if (info::ecf->input.omit_midpoints) {
		size_t shrinked = 0;
		for (size_t i = 0, offset = 0; i < etype.size(); ++i) {
			if (Mesh::edata[etype[i]].nodes != Mesh::edata[etype[i]].coarseNodes) {
				for (esint n = 0; n < Mesh::edata[etype[i]].coarseNodes; n++) {
					enodes[shrinked + n] = enodes[offset + n];
				}
				offset += Mesh::edata[etype[i]].nodes;
				shrinked += Mesh::edata[etype[i]].coarseNodes;
				esize[i] = Mesh::edata[etype[i]].coarseNodes;
			}
			switch (etype[i]) {
			case (int)Element::CODE::LINE3:     etype[i] = (int)Element::CODE::LINE2; break;
			case (int)Element::CODE::TRIANGLE6: etype[i] = (int)Element::CODE::TRIANGLE3; break;
			case (int)Element::CODE::SQUARE8:   etype[i] = (int)Element::CODE::SQUARE4; break;
			case (int)Element::CODE::TETRA10:   etype[i] = (int)Element::CODE::TETRA4; break;
			case (int)Element::CODE::PYRAMID13: etype[i] = (int)Element::CODE::PYRAMID5; break;
			case (int)Element::CODE::PRISMA15:  etype[i] = (int)Element::CODE::PRISMA6; break;
			case (int)Element::CODE::HEXA20:    etype[i] = (int)Element::CODE::HEXA8; break;
			default: break;
			}
		}
		enodes.resize(shrinked);
	}

	if (info::ecf->input.transformations.size()) {
		removeDuplicates = true;
		int instance = 1;
		Geometry geometry;
		selectGeometry(geometry);
		for (auto t = info::ecf->input.transformations.begin(); t != info::ecf->input.transformations.end(); ++t) {
			switch (t->second.transformation) {
			case InputTransformationConfiguration::TRANSFORMATION::TRANSLATE: translate(t->second, geometry, instance); break;
			case InputTransformationConfiguration::TRANSFORMATION::ROTATE: rotate(t->second, geometry, instance); break;
			case InputTransformationConfiguration::TRANSFORMATION::SCALE: scale(t->second, geometry, instance); break;
			case InputTransformationConfiguration::TRANSFORMATION::SHEAR: shear(t->second, geometry, instance); break;
			}
		}
	}
	profiler::synccheckpoint("transformation");

	for (auto nreg = info::ecf->input.node_regions.begin(); nreg != info::ecf->input.node_regions.end(); ++nreg) {
		auto &nodes = nregions[nreg->first];
		Point origin(nreg->second.x, nreg->second.y, nreg->second.z);
		Point end(nreg->second.x + nreg->second.lx, nreg->second.y + nreg->second.ly, nreg->second.z + nreg->second.lz);
		switch (nreg->second.shape) {
		case InputNodeRegionConfiguration::SHAPE::CIRCLE:
			for (size_t i = 0; i < coordinates.size(); ++i) {
				if ((origin - coordinates[i]).length() <= nreg->second.radius) {
					nodes.push_back(nIDs[i]);
				}
			}
			break;
		case InputNodeRegionConfiguration::SHAPE::BLOCK:
			for (size_t i = 0; i < coordinates.size(); ++i) {
				if (
						origin.x <= coordinates[i].x && coordinates[i].x <= end.x &&
						origin.y <= coordinates[i].y && coordinates[i].y <= end.y &&
						origin.z <= coordinates[i].z && coordinates[i].z <= end.z) {
					nodes.push_back(nIDs[i]);
				}
			}
			break;
		}
	}
	profiler::synccheckpoint("node_regions");

    for (auto noise = info::ecf->input.noise.cbegin(); noise != info::ecf->input.noise.end(); ++noise) {
        auto region = nregions.find(noise->first);
        if (region == nregions.end()) {
            eslog::globalerror("cannot set noise to the region '%s'. The region does not exist.\n", noise->first.c_str());
        } else {
            double &xx = noise->second.x.evaluator->getX();
            double &xy = noise->second.x.evaluator->getY();
            double &xz = noise->second.x.evaluator->getZ();
            double &yx = noise->second.x.evaluator->getX();
            double &yy = noise->second.x.evaluator->getY();
            double &yz = noise->second.x.evaluator->getZ();
            double &zx = noise->second.x.evaluator->getX();
            double &zy = noise->second.x.evaluator->getY();
            double &zz = noise->second.x.evaluator->getZ();
            for (size_t n = 0; n < region->second.size(); ++n) {
                const Point &p = coordinates[region->second[n]];
                xx = yx = zx = p.x;
                xy = yy = zy = p.y;
                xz = yz = zz = p.z;
                // to get deterministic noise to all processes
                std::srand(nIDs[region->second[n]]);
                int scale = std::rand() % 201 - 100;
                Point pp(noise->second.x.evaluator->evaluate(), noise->second.y.evaluator->evaluate(), noise->second.z.evaluator->evaluate());
                printf("%.3f %.3f %.3f >> ", pp.x, pp.y, pp.z);
                pp *= (scale / 100.);
                printf("%.3f %.3f %.3f (%d)\n", pp.x, pp.y, pp.z, scale);
                coordinates[region->second[n]] += Point(noise->second.x.evaluator->evaluate(), noise->second.y.evaluator->evaluate(), noise->second.z.evaluator->evaluate()) * (scale / 100.);
            }
        }

    }

    profiler::synccheckpoint("noise");

//	Communication::serialize([&] () {
////		{ std::ofstream os("XnIDs", std::ofstream::app); os << nIDs; os.close(); }
////		{ std::ofstream os("Xcoordinates", std::ofstream::app); os << coordinates; os.close(); }
////		{ std::ofstream os("XeIDs", std::ofstream::app); os << eIDs; os.close(); }
////		{ std::ofstream os("Xesize", std::ofstream::app); os << esize; os.close(); }
////		{ std::ofstream os("Xetype", std::ofstream::app); os << etype; os.close(); }
////		{ std::ofstream os("Xenodes", std::ofstream::app); os << enodes; os.close(); }
//		std::cout << nIDs.size() << " nIDs: " << nIDs;
//		std::cout << coordinates.size() << " coor: " << coordinates;
//
//		std::cout << eIDs.size() << " eIDs: " << eIDs;
//		std::cout << etype.size() << " etypes: " << etype;
//		std::cout << esize.size() << " esizes: " << esize;
//		std::cout << enodes.size() << " enodes: " << enodes;
//		std::cout << body.size() << " bodies: " << body;
//		std::cout << material.size() << " mats: " << material;
//
//		for (auto r = nregions.begin(); r != nregions.end(); ++r) {
//			std::cout << r->first << "[n]: " << r->second;
//		}
//		for (auto r = eregions.begin(); r != eregions.end(); ++r) {
//			std::cout << r->first << "[e]: " << r->second;
//		}
//	});

	switch (type) {
	case TYPE::GENERAL:
	case TYPE::SORTED: // TODO: optimized builder
		if (info::mpi::size > 1) {
			ScatteredInput{*this};
		} else {
			SequentialInput{*this};
		}
		break;
	case TYPE::GENERATED:
		GeneratedInput{*this, false};
	}
	info::mesh->orientation = orientation;

	const auto &selection = info::ecf->input.selection;
	for (auto sphere = selection.sphere.cbegin(); sphere != selection.sphere.cend(); ++sphere) {
	    Point center(sphere->second.cx, sphere->second.cy, sphere->second.cz);
	    std::vector<esint> nodes;
	    for (esint n = 0; n < info::mesh->nodes->size; ++n) {
	        if ((info::mesh->nodes->coordinates->datatarray()[n] - center).length() < sphere->second.r) {
	            nodes.push_back(n);
	        }
	    }
	    info::mesh->boundaryRegions.push_back(new BoundaryRegionStore(sphere->first));
	    info::mesh->boundaryRegions.back()->dimension = 1;
	    info::mesh->boundaryRegions.back()->nodes = new serializededata<esint, esint>(1, { (size_t)info::env::threads, nodes });
	}

	profiler::syncend("build_mesh");
}

void MeshBuilder::selectGeometry(Geometry &geometry)
{
	// TODO: generalize for an arbitrary selection
	std::vector<esint> max = { 0, 0 }, ids(2);
	if (nIDs.size()) {
		max[0] = *std::max_element(nIDs.begin(), nIDs.end());
	}
	if (eIDs.size()) {
		max[1] = *std::max_element(eIDs.begin(), eIDs.end());
	}
	Communication::allReduce(max.data(), ids.data(), 2, MPITools::getType<esint>().mpitype, MPI_MAX);
	geometry.nids = ids[0] + 1;
	geometry.eids = ids[1] + 1;

	geometry.nodes.first = 0;
	geometry.nodes.second = nIDs.size();
	geometry.elements.first = 0;
	geometry.elements.second = eIDs.size();
	geometry.enodes.first = 0;
	geometry.enodes.second = enodes.size();

	auto reg = nregions.begin();
	geometry.nregsize.resize(nregions.size());
	for (size_t r = 0; r < nregions.size(); ++r, ++reg) {
		geometry.nregsize[r].first = 0;
		geometry.nregsize[r].second = reg->second.size();
	}
	reg = eregions.begin();
	geometry.eregsize.resize(eregions.size());
	for (size_t r = 0; r < eregions.size(); ++r, ++reg) {
		geometry.eregsize[r].first = 0;
		geometry.eregsize[r].second = reg->second.size();
	}
}

template <typename TData>
static void _duplicate(std::vector<TData> &data, size_t begin, size_t end, esint offset)
{
	data.reserve(data.size() + (end - begin));
	for (size_t i = begin; i < end; ++i) {
		data.push_back(data[i] + offset);
	}
}

void MeshBuilder::duplicate(Geometry &source, int instance)
{
	_duplicate(nIDs, source.nodes.first, source.nodes.second, instance * source.nids);
	_duplicate(eIDs, source.elements.first, source.elements.second, instance * source.eids);
	_duplicate(etype, source.elements.first, source.elements.second, 0);
	_duplicate(esize, source.elements.first, source.elements.second, 0);
	_duplicate(enodes, source.enodes.first, source.enodes.second, instance * source.nids);
	_duplicate(body, source.elements.first, source.elements.second, 0);
	_duplicate(material, source.elements.first, source.elements.second, 0);

	if (_nranks.size()) {
		_nrankdist.reserve(_nrankdist.size() + (source.nodes.second - source.nodes.first));
		for (size_t i = source.nodes.first; i < source.nodes.second; ++i) {
			_nrankdist.push_back(_nrankdist.back() + _nrankdist[i + 1] - _nrankdist[i]);
			for (int r = _nrankdist[i]; r < _nrankdist[i + 1]; ++r) {
				_nranks.push_back(_nranks[r]);
			}
		}
	}

	auto reg = nregions.begin();
	for (size_t r = 0; r < nregions.size(); ++r, ++reg) {
		_duplicate(reg->second, source.nregsize[r].first, source.nregsize[r].second, instance * source.nids);
	}
	reg = eregions.begin();
	for (size_t r = 0; r < eregions.size(); ++r, ++reg) {
		_duplicate(reg->second, source.eregsize[r].first, source.eregsize[r].second, instance * source.eids);
	}
}

void MeshBuilder::translate(InputTransformationConfiguration &transformation, Geometry &source, int &instance)
{
	for (int i = 0; i < transformation.instances; ++i, ++instance) {
		duplicate(source, instance);
		coordinates.reserve(coordinates.size() + (source.nodes.second - source.nodes.first));
		for (size_t n = source.nodes.first; n < source.nodes.second; ++n) {
			coordinates.push_back(Point(
					coordinates[n].x + instance * transformation.x,
					coordinates[n].y + instance * transformation.y,
					coordinates[n].z + instance * transformation.z));
		}
	}
}

void MeshBuilder::rotate(InputTransformationConfiguration &transformation, Geometry &source, int &instance)
{
	for (int i = 0; i < transformation.instances; ++i, ++instance) {
		double cos = std::cos(instance * M_PI * transformation.z / 180), sin = std::sin(instance * M_PI * transformation.z / 180);
		duplicate(source, instance);
		coordinates.reserve(coordinates.size() + (source.nodes.second - source.nodes.first));
		for (size_t n = source.nodes.first; n < source.nodes.second; ++n) {
			coordinates.push_back(Point(
					cos * coordinates[n].x - sin * coordinates[n].y,
					sin * coordinates[n].x + cos * coordinates[n].y,
					coordinates[n].z));
		}
	}
}

void MeshBuilder::scale(InputTransformationConfiguration &transformation, Geometry &source, int &instance)
{

}

void MeshBuilder::shear(InputTransformationConfiguration &transformation, Geometry &source, int &instance)
{

}








