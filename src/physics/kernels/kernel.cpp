
#include "kernel.h"
#include "solverdataprovider/provider.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "physics/system/builder/builder.h"
#include "wrappers/mpi/communication.h"

#include <cmath>
#include <algorithm>

using namespace espreso;

Kernel::InstanceFiller::InstanceFiller(esint nvectors)
: begin(0), end(0), interval(-1), invalid(0), DOFs(0), insertK(true), insertM(true), insertC(true), insertR(true), insertF(true),
  reduction(1), K(NULL), M(NULL), C(NULL), R(NULL), F(NULL), offset(NULL)
{
	Re.initVectors(nvectors);
	Fe.initVectors(nvectors);
}

void Kernel::ParameterInfo::forceResize()
{
	isconstant.clear();
	isconstant.resize(inputs.size(), false);
	smartResize();
}

void Kernel::ParameterInfo::smartResize()
{
	if (data != NULL) {
		return;
	}

	for (size_t i = 0; i < inputs.size(); ++i) {
		bool constant = isconstant[i];
		for (auto it = inputs[i].begin(); it != inputs[i].end(); ++it) {
			constant &= it->parameter->constantInterval(i / dimensions);
		}
		isconstant[i] = constant;
	}

	nodes.resize(info::mesh->elements->eintervals.size(), 0);
	gps.resize(info::mesh->elements->eintervals.size(), 0);

	std::vector<std::vector<esint> > dist(info::env::threads);

	esint sum = 0, index = 0;
	for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei, ++index) {
		nodes[index] = Mesh::edata[ei->code].nodes;
		gps[index] = Mesh::edata[ei->code].N->size();
		dist[0].push_back(sum);
		if (constantInterval(index)) {
			switch (range) {
			case Range::AGREGATED:
			case Range::PER_NODE:
			case Range::PER_GP:
				sum += dimensions;
				break;
			case Range::PER_NODE_GP:
				sum += dimensions * Mesh::edata[ei->code].nodes * Mesh::edata[ei->code].N->size();
				break;
			case Range::PER_GP_GP:
				sum += dimensions * Mesh::edata[ei->code].N->size() * Mesh::edata[ei->code].N->size();
				break;
			case Range::EACH_NODE:
				sum += dimensions * Mesh::edata[ei->code].nodes;
				break;
			case Range::EACH_GP:
				sum += dimensions * Mesh::edata[ei->code].N->size();
				break;
			case Range::EACH_NODE_GP:
				sum += dimensions * Mesh::edata[ei->code].nodes * Mesh::edata[ei->code].N->size();
				break;
			}
		} else {
			switch (range) {
			case Range::AGREGATED:
				sum += (ei->end - ei->begin) * dimensions;
				break;
			case Range::PER_NODE:
				sum += (ei->end - ei->begin) * dimensions * Mesh::edata[ei->code].nodes;
				break;
			case Range::PER_GP:
				sum += (ei->end - ei->begin) * dimensions * Mesh::edata[ei->code].N->size();
				break;
			case Range::PER_NODE_GP:
				sum += (ei->end - ei->begin) * dimensions * Mesh::edata[ei->code].nodes * Mesh::edata[ei->code].N->size();
				break;
			case Range::PER_GP_GP:
				sum += (ei->end - ei->begin) * dimensions * Mesh::edata[ei->code].N->size() * Mesh::edata[ei->code].N->size();
				break;
			case Range::EACH_NODE:
				sum += (ei->end - ei->begin) * dimensions * Mesh::edata[ei->code].nodes;
				break;
			case Range::EACH_GP:
				sum += (ei->end - ei->begin) * dimensions * Mesh::edata[ei->code].N->size();
				break;
			case Range::EACH_NODE_GP:
				sum += (ei->end - ei->begin) * dimensions * Mesh::edata[ei->code].nodes * Mesh::edata[ei->code].N->size();
				break;
			}
		}
	}
	dist[0].push_back(sum);
	data = new serializededata<esint, double>(dist, tarray<double>(info::env::threads, sum, true));
}

void Kernel::ParameterInfo::eval()
{
	for (int d = 0; d < dimensions; ++d) {
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			if (evaluator[dimensions * i + d]) {
				Evaluator::Params params;
				insertParameterValues(params, i, d);
				evaluator[dimensions * i + d]->evalVector((data->begin() + i)->size() / dimensions, dimensions, params, (data->begin() + i)->data() + d);
			} else {
				for (size_t j = 0; j < (data->begin() + i)->size(); j += dimensions) {
					(data->begin() + i)->data()[j + d] = ivalues[dimensions * i + d];
				}
			}
		}
	}
}

void Kernel::ParameterInfo::swapInput(ParameterInfo &toreplace, ParameterInfo &replacement)
{
	for (size_t i = 0; i < inputs.size(); ++i) {
		for (size_t j = 0; j < inputs[i].size(); ++j) {
			if (inputs[i][j].parameter == &toreplace) {
				inputs[i][j].parameter = &replacement;
			}
		}
	}
}

void Kernel::ParameterInfo::insertParameterValues(Evaluator::Params &params, esint interval, int dimension)
{
	for (auto it = inputs[dimensions * interval + dimension].begin(); it != inputs[dimensions * interval + dimension].end(); ++it) {
		if (it->parameter->isconstant[it->parameter->dimensions * interval + it->offset]) {
			params.general.push_back({ (it->parameter->data->begin() + interval)->data(), it->offset, 0 });
		} else {
			params.general.push_back({ (it->parameter->data->begin() + interval)->data(), it->offset, it->parameter->dimensions });
		}
	}
}

Kernel::ParameterInfo::ParameterInfo(size_t intervals, int dimensions, ParameterInfo::Range range, double initial)
: range(range), constness(Status::EMPTY), dimensions(dimensions), data(NULL)
{
	ivalues.resize(dimensions * intervals, initial);
	isset.resize(dimensions * intervals, false);
	isconstant.resize(dimensions * intervals, true);
	evaluator.resize(dimensions * intervals, NULL);

	inputs.resize(dimensions * intervals);
}

Kernel::ElementParameterInfo::ElementParameterInfo(int dimensions, ParameterInfo::Range range, double initial)
: ParameterInfo(info::mesh->elements->eintervals.size(), dimensions, range, initial)
{

}

Kernel::BoundaryParameterInfo::BoundaryParameterInfo(int dimensions, ParameterInfo::Range range)
: ParameterInfo(info::mesh->boundaryRegions.size(), dimensions, range)
{

}

Kernel::InputParameterIterator::InputParameterIterator(ParameterInfo &info, esint interval, int increment)
: inc(info.isconstant[interval] ? 0 : increment), data((info.data->begin() + interval)->data())
{

}

Kernel::OuputParameterIterator::OuputParameterIterator(ParameterInfo &info, esint interval, int increment)
: inc(info.isconstant[interval] ? 0 : increment), data((info.data->begin() + interval)->data())
{

}

void Kernel::move(ElementParameterInfo &from, ElementParameterInfo &to)
{
	memcpy(to.data->datatarray().data(), from.data->datatarray().data(), from.data->datatarray().size() * sizeof(double));
}

void Kernel::move(NodeData *from, ElementParameterInfo &to)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; t++) {
		auto i = to.data->datatarray().begin(t);
		for (auto n = info::mesh->elements->procNodes->datatarray().cbegin(t); n != info::mesh->elements->procNodes->datatarray().cend(t); ++n) {
			for (int d = 0; d < from->dimension; ++d, ++i) {
				*i = from->data[*n * from->dimension + d];
			}
		}
	}
}

void Kernel::move(ElementParameterInfo &from, NodeData *to)
{
	if (from.constness == ParameterInfo::Status::GLOBAL) {
		auto i = from.data->datatarray().cbegin();
		for (size_t n = 0; n < to->data.size() / to->dimension; n++) {
			for (int d = 0; d < to->dimension; d++) {
				to->data[n * to->dimension + d] = i[d];
			}
		}
		return;
	}

	std::fill(to->data.begin(), to->data.end(), 0);

	int index = 0;
	for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei, ++index) {
		const double *value = (from.data->begin() + index)->data();
		int valueinc = from.isconstant[index] ? 0 : to->dimension;
		for (auto n = (info::mesh->elements->procNodes->begin() + ei->begin)->begin(); n != (info::mesh->elements->procNodes->begin() + ei->end)->begin(); ++n, value += valueinc) {
			for (int d = 0; d < to->dimension; d++) {
				to->data[*n * to->dimension + d] += value[d];
			}
		}
	}

	auto &nelements = info::mesh->nodes->elements->boundarytarray();
	for (size_t i = 0; i < to->data.size() / to->dimension; i++) {
		for (int d = 0; d < to->dimension; d++) {
			to->data[i * to->dimension + d] /= nelements[i + 1] - nelements[i];
		}
	}

	std::vector<std::vector<double> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());

	auto nranks = info::mesh->nodes->ranks->begin();
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		if (nranks->size() > 1) {
			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbors[noffset] < *r) {
						++noffset;
					}
					for (int d = 0; d < to->dimension; d++) {
						sBuffer[noffset].push_back(to->data[n * to->dimension + d]);
					}
				}
			}
		}
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
		rBuffer[n].resize(sBuffer[n].size());
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::error("ESPRESO internal error: exchange diagonal values.\n");
	}

	nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> nindex(info::mesh->neighbors.size());
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		if (nranks->size() > 1) {
			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbors[noffset] < *r) {
						++noffset;
					}
					for (int d = 0; d < to->dimension; d++) {
						to->data[n * to->dimension + d] += rBuffer[noffset][nindex[noffset]++];
					}
				}
			}
		}
	}
}

void Kernel::move(BoundaryParameterInfo &info, NodeData *data)
{
	// TODO: average crossing values
	int index = 0;
	for (size_t i = 0; i < info.isset.size(); ++i) {
		if (info.isset[i]) {
			const BoundaryRegionStore *bregion = info::mesh->boundaryRegions[index];
			const double *value = (info.data->begin() + index)->data();
			int valueinc = info.isconstant[index] ? 0 : data->dimension;
			for (auto n = bregion->nodes->datatarray().begin(); n != bregion->nodes->datatarray().end(); ++n, value += valueinc) {
				for (int d = 0; d < data->dimension; d++) {
					data->data[*n * data->dimension + d] += value[d];
				}
			}
		}
	}
}

void Kernel::smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order) const
{
	value = std::max(0.0, std::min((value - edge0) / (edge1 - edge0), 1.0));

	switch (order) {
	case 0:
		smoothStep = value;
		if (value == 0 || value == 1) {
			derivation = 0;
		} else {
			derivation = 1;
		}
		break;
	case 1:
		smoothStep = -2 * pow(value, 3) + 3 * pow(value, 2);
		derivation = -6 * pow(value, 2) + 6 * value;
		break;
	case 2:
		smoothStep = 6 * pow(value, 5) - 15 * pow(value, 4) + 10 * pow(value, 3);
		derivation = 30 * pow(value, 4) - 60*  pow(value, 3) + 30 * pow(value, 2);
		break;
	case 3:
		smoothStep = -20 * pow(value, 7) + 70 * pow(value, 6) - 84 * pow(value, 5) + 35 * pow(value, 4);
		derivation = -140 * pow(value, 6) + 420 * pow(value, 5) - 420 * pow(value, 4) + 140 * pow(value, 3);
		break;
	case 4:
		smoothStep = 70 * pow(value, 9) - 315 * pow(value, 8) + 540 * pow(value, 7) - 420 * pow(value, 6) + 126 * pow(value, 5);
		derivation = 630 * pow(value, 8) - 2520 * pow(value, 7) + 3780 * pow(value, 6) - 2520 * pow(value, 5) + 630 * pow(value, 4);
		break;
	case 5:
		smoothStep = -252 * pow(value, 11) + 1386 * pow(value, 10) - 3080 * pow(value, 9) + 3465 * pow(value, 8) - 1980 * pow(value, 7) + 462 * pow(value, 6);
		derivation = -2772 * pow(value, 10) + 13860 * pow(value, 9) - 27720 * pow(value, 8) + 27720 * pow(value, 7) - 13860 * pow(value, 6) + 2772 * pow(value, 5);
		break;
	case 6:
		smoothStep = 924 * pow(value, 13) - 6006 * pow(value, 12) + 16380 * pow(value, 11) - 24024 * pow(value, 10) + 20020 * pow(value, 9) - 9009 * pow(value, 8) + 1716 * pow(value, 7);
		derivation = 12012 * pow(value, 12) - 72072 * pow(value, 11) + 180180 * pow(value, 10) - 240240 * pow(value, 9) + 180180 * pow(value, 8) - 72072 * pow(value, 7) + 12012 * pow(value, 6);
		break;
	}

	derivation /= edge1 - edge0;
}


