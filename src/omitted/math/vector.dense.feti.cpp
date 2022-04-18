
#include "vector.dense.feti.h"
#include "vector.dense.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"

#include <cstddef>
#include <limits>

using namespace espreso;

VectorDenseFETI::VectorDenseFETI()
{

}

VectorDenseFETI::~VectorDenseFETI()
{

}

VectorDenseFETI* VectorDenseFETI::copy()
{
	return new VectorDenseFETI();
}

Vector* VectorDenseFETI::create()
{
	return new VectorDense();
}

double VectorDenseFETI::norm()
{
	if (this->duplications == DataDecomposition::DUPLICATION::SPLIT) {
		allGather();
		double square = 0;
		esint i = 0;
		for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
			double value = 0;
			if (ismy(map->begin()->domain)) {
				if (map->size() > 1) {
					for (size_t j = i; j < i + map->size(); ++j) {
						value += gathered[j];
					}
				} else {
					value = at(map->begin()->domain - doffset)->vals[map->begin()->index];
				}
				square += value * value;
			}
			if (map->size() > 1) {
				i += map->size();
			}
		}

		double norm;
		Communication::allReduce(&square, &norm, 1, MPI_DOUBLE, MPI_SUM);
		return std::sqrt(norm);
	} else {
		double square = 0;
		for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
			if (ismy(map->begin()->domain)) {
				square += at(map->begin()->domain - doffset)->vals[map->begin()->index] * at(map->begin()->domain - doffset)->vals[map->begin()->index];
			}
		}
		double norm;
		Communication::allReduce(&square, &norm, 1, MPI_DOUBLE, MPI_SUM);
		return std::sqrt(norm);
	}
}

double VectorDenseFETI::max()
{
	if (this->duplications == DataDecomposition::DUPLICATION::SPLIT) {
		eslog::internalFailure("implement max of SPLIT FETI vectors.\n");
	}

	double max = std::numeric_limits<double>::min();
	for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
		if (ismy(map->begin()->domain)) {
			if (max < at(map->begin()->domain - doffset)->vals[map->begin()->index]) {
				max = at(map->begin()->domain - doffset)->vals[map->begin()->index];
			}
		}
	}
	double gmax;
	Communication::allReduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX);
	return gmax;
}

double VectorDenseFETI::absmax()
{
	if (this->duplications == DataDecomposition::DUPLICATION::SPLIT) {
		eslog::internalFailure("implement max of SPLIT FETI vectors.\n");
	}

	double max = 0;
	for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
		if (ismy(map->begin()->domain)) {
			if (max < std::fabs(at(map->begin()->domain - doffset)->vals[map->begin()->index])) {
				max = std::fabs(at(map->begin()->domain - doffset)->vals[map->begin()->index]);
			}
		}
	}
	double gmax;
	Communication::allReduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX);
	return gmax;
}

double VectorDenseFETI::dot(const Vector *other)
{
	double res = 0;
	const VectorDenseFETI *_other = other->downcast<VectorDenseFETI>();

	if (
			this->duplications == DataDecomposition::DUPLICATION::SPLIT ||
			_other->duplications == DataDecomposition::DUPLICATION::SPLIT) {

		this->allGather();
		_other->allGather();

		double dot = 0;
		esint i = 0;
		for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
			if (ismy(map->begin()->domain)) {
				double a = 0;
				double b = 0;
				if (map->size() > 1) {
					for (size_t j = i; j < i + map->size(); ++j) {
						a += gathered[j];
						b += _other->gathered[j];
					}
				} else {
					a = at(map->begin()->domain - doffset)->vals[map->begin()->index];
					b = _other->at(map->begin()->domain - doffset)->vals[map->begin()->index];
				}
				dot += a * b;
			}
			if (map->size() > 1) {
				i += map->size();
			}
		}

		Communication::allReduce(&dot, &res, 1, MPI_DOUBLE, MPI_SUM);
	}

	if (
			this->duplications == DataDecomposition::DUPLICATION::DUPLICATE ||
			_other->duplications == DataDecomposition::DUPLICATION::SPLIT) {

		_other->allGather();

		double dot = 0;
		esint i = 0;
		for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
			if (ismy(map->begin()->domain)) {
				double a = at(map->begin()->domain - doffset)->vals[map->begin()->index];
				double b = 0;
				if (map->size() > 1) {
					for (size_t j = i; j < i + map->size(); ++j) {
						b += _other->gathered[j];
					}
				} else {
					b = _other->at(map->begin()->domain - doffset)->vals[map->begin()->index];
				}
				dot += a * b;
			}
			if (map->size() > 1) {
				i += map->size();
			}
		}

		Communication::allReduce(&dot, &res, 1, MPI_DOUBLE, MPI_SUM);
	}

	if (
			this->duplications == DataDecomposition::DUPLICATION::SPLIT ||
			_other->duplications == DataDecomposition::DUPLICATION::DUPLICATE) {

		this->allGather();

		double dot = 0;
		esint i = 0;
		for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
			if (ismy(map->begin()->domain)) {
				double a = 0;
				double b = _other->at(map->begin()->domain - doffset)->vals[map->begin()->index];
				if (map->size() > 1) {
					for (size_t j = i; j < i + map->size(); ++j) {
						a += gathered[j];
					}
				} else {
					a = at(map->begin()->domain - doffset)->vals[map->begin()->index];
				}
				dot += a * b;
			}
			if (map->size() > 1) {
				i += map->size();
			}
		}

		Communication::allReduce(&dot, &res, 1, MPI_DOUBLE, MPI_SUM);
	}

	if (
			this->duplications == DataDecomposition::DUPLICATION::DUPLICATE &&
			_other->duplications == DataDecomposition::DUPLICATION::DUPLICATE) {

		double dot = 0;
		for (auto map = dmap->cbegin(); map != dmap->cend(); ++map) {
			if (ismy(map->begin()->domain)) {
				dot += at(map->begin()->domain - doffset)->vals[map->begin()->index] * _other->at(map->begin()->domain - doffset)->vals[map->begin()->index];
			}
		}

		Communication::allReduce(&dot, &res, 1, MPI_DOUBLE, MPI_SUM);
	}
	return res;
}

void VectorDenseFETI::allGather() const
{
	DataDecomposition::allGather(*this);
}

void VectorDenseFETI::averageDuplications()
{
	allGather();
	auto map = dmap->begin();
	for (esint n = 0, i = 0, prev = 0; n < nshared; prev = shared[n++]) {
		map += shared[n] - prev;
		double sum = 0;
		for (size_t j = i; j < i + map->size(); ++j) {
			sum += gathered[j];
		}
		for (auto di = map->begin(); di != map->end(); ++di, ++i) {
			if (ismy(di->domain)) {
				at(di->domain - doffset)->vals[di->index] = sum / map->size();
			}
		}
	}
}

void VectorDenseFETI::sumDuplications()
{
	allGather();
	auto map = dmap->begin();
	for (esint n = 0, i = 0, prev = 0; n < nshared; prev = shared[n++]) {
		map += shared[n] - prev;
		double sum = 0;
		for (size_t j = i; j < i + map->size(); ++j) {
			sum += gathered[j];
		}
		for (auto di = map->begin(); di != map->end(); ++di, ++i) {
			if (ismy(di->domain)) {
				at(di->domain - doffset)->vals[di->index] = sum;
			}
		}
	}
}

VectorsDenseFETI::VectorsDenseFETI()
{
	initVectors(0);
}


VectorsDenseFETI::~VectorsDenseFETI()
{

}

VectorsDenseFETI* VectorsDenseFETI::copy()
{
	return new VectorsDenseFETI();
}

Vector* VectorsDenseFETI::create()
{
	return new VectorDenseFETI();
}

void VectorsDenseFETI::initDomains(DataDecomposition::DUPLICATION duplications, esint domains)
{
	holder()->initDomains(duplications, domains);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->initDomains(duplications, domains);
	}
}

void VectorsDenseFETI::resizeDomain(esint domain, esint size)
{
	holder()->at(domain)->resize(size * nvectors);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->at(domain)->shallowCopyFromHolder(holder()->at(domain), n, nvectors);
	}
}

void VectorsDenseFETI::fillDecomposition(esint rank, esint ranks, esint nneighbors, esint *distribution, int *neighbors, const serializededata<esint, DI> *dmap)
{
	holder()->fillDecomposition(rank, ranks, nneighbors, distribution, neighbors, dmap);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->DataDecomposition::shallowCopy(holder());
	}
}

void VectorsDenseFETI::averageDuplications()
{
	// TODO: implement more effective way (through dmap)
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->averageDuplications();
	}
}

void VectorsDenseFETI::sumDuplications()
{
	// TODO: implement more effective way (through dmap)
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->sumDuplications();
	}
}



