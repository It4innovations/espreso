
#ifndef SRC_INPUT_INPUT_H_
#define SRC_INPUT_INPUT_H_

#include "basis/containers/allocators.h"
#include "basis/containers/point.h"
#include "mesh/element.h"

#include <string>
#include <vector>

namespace espreso {

struct InputConfiguration;
class Mesh;

class Input {
public:
	virtual void load(const InputConfiguration &configuration) = 0;
	virtual void build(Mesh &mesh) = 0;

	virtual ~Input() {}
};


/*
 * Element connectivity is described from offsets to coordinates.
 * Coordinates are continuous.
 */

struct OrderedMeshDatabase {

	struct Offset {
		esint offset, start, size;
	};

	struct Region {
		std::string name;
		esint offset, size; // begin, end offset to region
	};

	struct Values {
		struct ValuesArray {
			ValuesArray(esint begin, int dimension): offset(begin), dimension(dimension) {}

			esint offset;
			int dimension; // -1 for all dimensions
			std::vector<esfloat> values;
		};

		Values(const std::string &name, int dimension): name(name), dimension(dimension) {}

		std::string name;
		int dimension;
		std::vector<ValuesArray> values;
	};

	std::vector<_Point<esfloat>, initless_allocator<_Point<esfloat> > > coordinates;

	std::vector<Element::CODE, initless_allocator<Element::CODE> > etype;
	std::vector<esint, initless_allocator<esint> > enodes;

	std::vector<Offset> noffsets, eoffsets;
	std::vector<Region> nregions, eregions;
	std::vector<Values> nvalues, evalues;

	static bool chunk(const esint &mpichunk, const int &rank, const std::vector<Offset> &offsets, std::vector<Offset>::const_iterator &it, esint &begin, esint &end)
	{
		if (it != offsets.end() && end == it->offset + it->size) {
			++it;
		}
		if (it == offsets.end()) {
			return false;
		}
		if (it->offset + end - it->start == mpichunk * (rank + 1)) {
			return false;
		}
		begin = end = 0;
		while (it != offsets.end() && it->offset + it->size < mpichunk * rank) { ++it; }
		if (it != offsets.end() && it->offset < mpichunk * (rank + 1)) {
			begin = std::max(it->offset, mpichunk * rank);
			end = std::min(it->offset + it->size, mpichunk * (rank + 1));
			begin = it->start + begin - it->offset;
			end = it->start + end - it->offset;
		}
		return begin != end;
	}

	void clearNodes()
	{
		std::vector<Offset>().swap(noffsets);
		std::vector<_Point<esfloat>, initless_allocator<_Point<esfloat> > >().swap(coordinates);
	}

	void clearElements()
	{
		std::vector<Offset>().swap(eoffsets);
		std::vector<Element::CODE, initless_allocator<Element::CODE> >().swap(etype);
		std::vector<esint, initless_allocator<esint> >().swap(enodes);
	}

	void clearRegions()
	{
		std::vector<Region>().swap(nregions);
		std::vector<Region>().swap(eregions);
	}

	void clearValues()
	{
		std::vector<Values>().swap(nvalues);
		std::vector<Values>().swap(evalues);
	}
};

inline bool operator<(const OrderedMeshDatabase::Offset &o1, const OrderedMeshDatabase::Offset &o2) { return o1.offset < o2.offset; }

}

#endif /* SRC_INPUT_INPUT_H_ */
