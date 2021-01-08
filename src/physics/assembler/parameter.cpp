
#include "parameter.h"

#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "math/matrix.dense.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

ParameterData::ParameterData(PerElementSize mask, int intervals): size(mask), data(NULL), isset(false)
{
	isconst.resize(intervals, true);
	version.resize(intervals, 0);
}

ElementParameterData::ElementParameterData(PerElementSize mask)
: ParameterData(mask, intervals())
{

}

int ElementParameterData::intervals()
{
	return info::mesh->elements->eintervals.size();
}

BoundaryParameterData::BoundaryParameterData(int region, PerElementSize mask)
: ParameterData(mask, intervals(region)), region(region)
{

}

int BoundaryParameterData::regions()
{
	return info::mesh->boundaryRegions.size();
}

int BoundaryParameterData::intervals(int region)
{
	return std::max(1, (int)info::mesh->boundaryRegions[region]->eintervals.size());
}

BoundaryParameterPack::BoundaryParameterPack(PerElementSize mask)
: _mask(mask)
{
	regions.reserve(info::mesh->boundaryRegions.size());
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		regions.emplace_back(BoundaryParameterData(r, _mask));
	}
}

void ParameterData::addInput(const ECFExpression &exp, int interval)
{

}

void ParameterData::addInput(const ParameterData &p)
{
	for (size_t i = 0; i < isconst.size(); ++i) {
		isconst[i] = isconst[i] && p.isconst[i];
	}
}

// currently we assume that serialized data are never updated
void ParameterData::addInput(const serializededata<esint, esint>* p)
{
	setConstness(false);
}

void ParameterData::addInput(const serializededata<esint, Point>* p)
{
	setConstness(false);
}

void ParameterData::addInput(const NodeData* p)
{
	setConstness(false);
}

void ParameterData::addInput(int interval, const serializededata<esint, Point>* p)
{
	isconst[interval] = false;
}

void ParameterData::setConstness(bool constness)
{
	std::fill(isconst.begin(), isconst.end(), constness);
}

void ElementParameterData::resize()
{
	std::vector<std::vector<esint> > distribution(info::env::threads);

	distribution[0].push_back(0);
	for (int t = 0; t < info::env::threads; ++t) {
		esint sum = 0;
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				esint isize = increment(i);
				if (!isconst[i]) {
					isize *= info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
				}
				distribution[t].push_back(isize + sum);
				sum += isize;
			}
		}
	}
	utils::threadDistributionToFullDistribution(distribution);

	std::vector<size_t> datadistribution;
	for (int t = 0; t < info::env::threads; ++t) {
		datadistribution.push_back(distribution[t].front());
	}
	datadistribution.push_back(distribution.back().back());

	if (data) {
		delete data;
	}
	data = new serializededata<esint, double>(distribution, tarray<double>(datadistribution, 1));
	isset = true;
}

int ElementParameterData::increment(int interval) const
{
	return
			size.n *
			std::pow(info::mesh->dimension, size.ndimension) *
			std::pow(info::mesh->edata[info::mesh->elements->eintervals[interval].code].dimension, size.edimension) *
			std::pow(info::mesh->edata[info::mesh->elements->eintervals[interval].code].nodes, size.node) *
			std::pow(info::mesh->edata[info::mesh->elements->eintervals[interval].code].N->size(), size.gp);
}

int ElementParameterData::increment(PerElementSize size, int interval) const
{
	return
			size.n *
			std::pow(info::mesh->dimension, size.ndimension) *
			std::pow(info::mesh->edata[info::mesh->elements->eintervals[interval].code].dimension, size.edimension) *
			std::pow(info::mesh->edata[info::mesh->elements->eintervals[interval].code].nodes, size.node) *
			std::pow(info::mesh->edata[info::mesh->elements->eintervals[interval].code].N->size(), size.gp);
}

void BoundaryParameterData::resize()
{
	if (data) {
		delete data;
	}
	if (region == 0) {
		// ALL_NODES are kept empty (is it possible to have boundary condition on all nodes?)
		data = new serializededata<esint, double>({ 0, 0 }, tarray<double>(0, 0));
		return;
	}
	if (info::mesh->boundaryRegions[region]->dimension == 0) {
		esint dimension = size.n * std::pow(info::mesh->dimension, size.ndimension);
		if (isconst[0]) {
			data = new serializededata<esint, double>(dimension, tarray<double>(1, dimension));
		} else {
			data = new serializededata<esint, double>(dimension, tarray<double>(info::mesh->boundaryRegions[region]->nodes->datatarray().distribution(), dimension));
		}
		return;
	}
	std::vector<std::vector<esint> > distribution(info::env::threads);

	distribution[0].push_back(0);
	for (int t = 0; t < info::env::threads; ++t) {
		esint sum = 0;
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			for (esint i = info::mesh->boundaryRegions[region]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[region]->eintervalsDistribution[d + 1]; ++i) {
				esint isize = increment(i);
				if (!isconst[i]) {
					isize *= info::mesh->boundaryRegions[region]->eintervals[i].end - info::mesh->boundaryRegions[region]->eintervals[i].begin;
				}
				distribution[t].push_back(isize + sum);
				sum += isize;
			}
		}
	}
	utils::threadDistributionToFullDistribution(distribution);

	std::vector<size_t> datadistribution;
	for (int t = 0; t < info::env::threads; ++t) {
		datadistribution.push_back(distribution[t].front());
	}
	datadistribution.push_back(distribution.back().back());

	data = new serializededata<esint, double>(distribution, tarray<double>(datadistribution, 1));
}

int BoundaryParameterData::increment(int interval) const
{
	return
			size.n *
			std::pow(info::mesh->dimension, size.ndimension) *
			std::pow(info::mesh->edata[info::mesh->boundaryRegions[region]->eintervals[interval].code].dimension, size.edimension) *
			std::pow(info::mesh->edata[info::mesh->boundaryRegions[region]->eintervals[interval].code].nodes, size.node) *
			std::pow(info::mesh->edata[info::mesh->boundaryRegions[region]->eintervals[interval].code].N->size(), size.gp);
}

int BoundaryParameterData::increment(PerElementSize size, int interval) const
{
	return
			size.n *
			std::pow(info::mesh->dimension, size.ndimension) *
			std::pow(info::mesh->edata[info::mesh->boundaryRegions[region]->eintervals[interval].code].dimension, size.edimension) *
			std::pow(info::mesh->edata[info::mesh->boundaryRegions[region]->eintervals[interval].code].nodes, size.node) *
			std::pow(info::mesh->edata[info::mesh->boundaryRegions[region]->eintervals[interval].code].N->size(), size.gp);
}
