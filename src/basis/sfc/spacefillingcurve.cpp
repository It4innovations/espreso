
#include "spacefillingcurve.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"

#include <utility>
#include <limits>
#include <algorithm>

namespace espreso {

template <typename T>
static void _build(SpaceFillingCurve<T> *sfc, size_t dimension, size_t depth, size_t npoints, _Point<T>* coordinates)
{
	if (dimension != 2 && dimension != 3) {
		eslog::globalerror("incorrect mesh dimension ='%ld'.\n", dimension);
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	double dmax = std::numeric_limits<double>::max();
	std::vector<_Point<T> > mins(threads, Point(dmax, dmax, dmax)), maxs(threads, _Point<T>(-dmax, -dmax, -dmax));

	std::vector<size_t> cdistribution = tarray<size_t>::distribute(threads, npoints);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		_Point<T> tmin = _Point<T>(dmax, dmax, dmax), tmax = _Point<T>(-dmax, -dmax, -dmax);
		for (size_t i = cdistribution[t]; i < cdistribution[t + 1]; i++) {
			tmin.x = std::min(tmin.x, coordinates[i].x);
			tmin.y = std::min(tmin.y, coordinates[i].y);
			tmin.z = std::min(tmin.z, coordinates[i].z);
			tmax.x = std::max(tmax.x, coordinates[i].x);
			tmax.y = std::max(tmax.y, coordinates[i].y);
			tmax.z = std::max(tmax.z, coordinates[i].z);
		}

		mins[t] = tmin;
		maxs[t] = tmax;
	}

	for (size_t t = 1; t < threads; t++) {
		mins[0].x = std::min(mins[t].x, mins[0].x);
		mins[0].y = std::min(mins[t].y, mins[0].y);
		mins[0].z = std::min(mins[t].z, mins[0].z);
		maxs[0].x = std::max(maxs[t].x, maxs[0].x);
		maxs[0].y = std::max(maxs[t].y, maxs[0].y);
		maxs[0].z = std::max(maxs[t].z, maxs[0].z);
	}

	Communication::allReduce(&mins[0].x, &sfc->origin.x, 3, MPITools::getType<esfloat>().mpitype, MPI_MIN);
	Communication::allReduce(&maxs[0].x, &sfc->size.x, 3, MPITools::getType<esfloat>().mpitype, MPI_MAX);

	sfc->size -= sfc->origin - Point(1e-6, 1e-6, 1e-6);
}

template <typename T>
static void _addSFCNeighbors(const SpaceFillingCurve<T> *sfc, size_t depth, size_t index, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors)
{
	size_t x, y, z = 0, nsize;
	if (sfc->dimension == 2) {
		sfc->D1toD2((size_t)1 << depth, index, x, y);
		nsize = neighbors.size();
		sfc->addXYNeighbors(depth, x, y, splitters, neighbors);
		for (size_t i = nsize; i < neighbors.size(); i++) {
			size_t n = (size_t)1 << neighbors[i].first;
			neighbors[i].second = sfc->D2toD1(n, neighbors[i].second % n, neighbors[i].second / n);
		}
	}
	if (sfc->dimension == 3) {
		sfc->D1toD3((size_t)1 << depth, index, x, y, z);
		nsize = neighbors.size();
		sfc->addXYZNeighbors(depth, x, y, z, splitters, neighbors);
		for (size_t i = nsize; i < neighbors.size(); i++) {
			size_t n = (size_t)1 << neighbors[i].first;
			neighbors[i].second = sfc->D3toD1(n, neighbors[i].second % n, neighbors[i].second % (n * n) / n, neighbors[i].second / (n * n));
		}
	}
}

static std::pair<size_t, size_t> _getXYZBucket(size_t depth, size_t x, size_t y, size_t z)
{
	std::pair<size_t, size_t> bucket;

	size_t coarsenig = (size_t)1 << depth;
	size_t cx, cy, cz;

	size_t d = 0, n = 1;
	do {
		cx = x / coarsenig;
		cy = y / coarsenig;
		cz = z / coarsenig;

		bucket.first = d;
		bucket.second = cz * n * n + cy * n + cx;

		coarsenig /= 2; n *= 2;
//	} while (++d <= depth && std::binary_search(_refinedxyz[bucket.first].begin(), _refinedxyz[bucket.first].end(), bucket.second));
	} while (++d <= depth);

	return bucket;
}

template <typename T>
static void _addXYNeighbors(const SpaceFillingCurve<T> *sfc, size_t depth, size_t x, size_t y, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors)
{
	std::vector<std::pair<size_t, size_t> > potential;

	size_t n = (size_t)1 << depth;
	size_t nn, xx, yy;

	for (int ox = -1; ox <= 1; ox++) {
		for (int oy = -1; oy <= 1; oy++) {
			if (x + ox < n && y + oy < n) {

				potential.clear();
				potential.push_back(_getXYZBucket(depth, x + ox, y + oy, 0));

				// neighbors should be finer
				for (size_t i = 0; i < potential.size(); i++) {
					size_t n = (size_t)1 << potential[i].first;
					size_t bucket = sfc->D2toD1(n, potential[i].second % n, potential[i].second / n);
					size_t bstep = sfc->buckets(sfc->depth) / sfc->buckets(potential[i].first);
					esint first = bucket * bstep;
					esint last  = bucket * bstep + bstep;

					auto sit = std::lower_bound(splitters.begin(), splitters.end(), first + 1);
					if (*sit < last) {
						nn = (size_t)1 << potential[i].first;
						xx = potential[i].second % nn;
						yy = potential[i].second / nn;
						nn = nn << 1;
						for (size_t oox = (ox == -1 ? 1 : 0); oox < (ox == 1 ? 1 : 2); oox++) {
							for (size_t ooy = (oy == -1 ? 1 : 0); ooy < (oy == 1 ? 1 : 2); ooy++) {
								potential.push_back(std::pair<size_t, size_t>(potential[i].first + 1, nn * (2 * yy + ooy) + 2 * xx + oox));
							}
						}
					} else {
						neighbors.push_back(potential[i]);
					}
				}

			}
		}
	}
}

template <typename T>
static void _addXYZNeighbors(const SpaceFillingCurve<T> *sfc, size_t depth, size_t x, size_t y, size_t z, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors)
{
	std::vector<std::pair<size_t, size_t> > potential;

	size_t n = (size_t)1 << depth;
	size_t nn, xx, yy, zz;

	for (int ox = -1; ox <= 1; ox++) {
		for (int oy = -1; oy <= 1; oy++) {
			for (int oz = -1; oz <= 1; oz++) {
				if (x + ox < n && y + oy < n && z + oz < n) {

					potential.clear();
					potential.push_back(_getXYZBucket(depth, x + ox, y + oy, z + oz));

					for (size_t i = 0; i < potential.size(); i++) {
						size_t n = (size_t)1 << potential[i].first;
						size_t bucket = sfc->D3toD1(n, potential[i].second % n, potential[i].second % (n * n) / n, potential[i].second / (n * n));
						size_t bstep = sfc->buckets(sfc->depth) / sfc->buckets(potential[i].first);
						esint first = bucket * bstep;
						esint last  = bucket * bstep + bstep;

						auto sit = std::lower_bound(splitters.begin(), splitters.end(), first + 1); // +1 since there can be invalid check if we hit interval begin
						if (*sit < last) {
							nn = (size_t)1 << potential[i].first;
							xx = potential[i].second % nn;
							yy = potential[i].second % (nn * nn) / nn;
							zz = potential[i].second / (nn * nn);
							nn = nn << 1;
							for (size_t oox = (ox == -1 ? 1 : 0); oox < (ox == 1 ? 1 : 2); oox++) {
								for (size_t ooy = (oy == -1 ? 1 : 0); ooy < (oy == 1 ? 1 : 2); ooy++) {
									for (size_t ooz = (oz == -1 ? 1 : 0); ooz < (oz == 1 ? 1 : 2); ooz++) {
										potential.push_back(std::pair<size_t, size_t>(potential[i].first + 1, nn * nn * (2 * zz + ooz) + nn * (2 * yy + ooy) + 2 * xx + oox));
									}
								}
							}
						} else {
							neighbors.push_back(potential[i]);
						}
					}

				}
			}
		}
	}
}

template <>
SpaceFillingCurve<float>::SpaceFillingCurve(size_t dimension, size_t depth, size_t npoints, _Point<float>* coordinates)
: dimension(dimension), depth(depth), n((size_t)1 << depth)
{
	_build(this, dimension, depth, npoints, coordinates);
}

template <>
SpaceFillingCurve<double>::SpaceFillingCurve(size_t dimension, size_t depth, size_t npoints, _Point<double>* coordinates)
: dimension(dimension), depth(depth), n((size_t)1 << depth)
{
	_build(this, dimension, depth, npoints, coordinates);
}

template <> void SpaceFillingCurve<float>::addSFCNeighbors(size_t depth, size_t index, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const { _addSFCNeighbors(this, depth, index, splitters, neighbors); }
template <> void SpaceFillingCurve<double>::addSFCNeighbors(size_t depth, size_t index, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const { _addSFCNeighbors(this, depth, index, splitters, neighbors); }

template <> void SpaceFillingCurve<float>::addXYNeighbors(size_t depth, size_t x, size_t y, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const { _addXYNeighbors(this, depth, x, y, splitters, neighbors); }
template <> void SpaceFillingCurve<double>::addXYNeighbors(size_t depth, size_t x, size_t y, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const { _addXYNeighbors(this, depth, x, y, splitters, neighbors); }

template <> void SpaceFillingCurve<float>::addXYZNeighbors(size_t depth, size_t x, size_t y, size_t z, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const { _addXYZNeighbors(this, depth, x, y, z, splitters, neighbors); }
template <> void SpaceFillingCurve<double>::addXYZNeighbors(size_t depth, size_t x, size_t y, size_t z, const std::vector<esint> &splitters, std::vector<std::pair<size_t, size_t> > &neighbors) const { _addXYZNeighbors(this, depth, x, y, z, splitters, neighbors); }

}
