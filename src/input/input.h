
#ifndef SRC_INPUT_INPUT_H_
#define SRC_INPUT_INPUT_H_

#include "basis/containers/allocators.h"
#include "basis/containers/point.h"

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

	struct Region {
		std::string name;
		esint begin, end; // start, end offset to region
	};

	struct Values {
		struct ValuesArray {
			ValuesArray(esint begin, int dimension): begin(begin), dimension(dimension) {}

			esint begin;
			int dimension; // -1 for all dimensions
			std::vector<esfloat> values;
		};

		Values(const std::string &name, int dimension): name(name), dimension(dimension) {}

		std::string name;
		int dimension;
		std::vector<ValuesArray> values;
	};

	std::vector<esint> noffset;
	std::vector<_Point<esfloat> > coordinates;

	std::vector<esint> eoffset;
	std::vector<int> etype;
	std::vector<int> esize;
	std::vector<esint> enodes;

	std::vector<Region> nregions, eregions;
	std::vector<Values> nvalues, evalues;

	void clearNodes()
	{
		std::vector<esint>().swap(noffset);
		std::vector<_Point<esfloat> >().swap(coordinates);
	}

	void clearElements()
	{
		std::vector<esint>().swap(eoffset);
		std::vector<int>().swap(etype);
		std::vector<int>().swap(esize);
		std::vector<esint>().swap(enodes);
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

}

#endif /* SRC_INPUT_INPUT_H_ */
