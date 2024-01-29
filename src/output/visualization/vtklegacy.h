
#ifndef SRC_OUTPUT_VISUALIZATION_COLLECTED_VTKLEGACY_H_
#define SRC_OUTPUT_VISUALIZATION_COLLECTED_VTKLEGACY_H_

#include "basis/containers/allocators.h"
#include "visualization.h"
#include "writer/vtkwritter.h"
#include <string>

namespace espreso {

struct RegionStore;
struct ElementsRegionStore;
struct BoundaryRegionStore;

class VTKLegacy: public Visualization {
public:
	VTKLegacy();
	~VTKLegacy();

	void updateMesh();
	void updateMonitors(const step::Step &step);
	void updateSolution(const step::Step &step, const step::Time &time);
	void updateSolution(const step::Step &step, const step::Frequency &frequency);

protected:
	void updateSolution(const step::Step &step, const std::string &dir, const std::string &name, const std::string &suffix);

	void insertHeader();
	void insertPoints(const RegionStore *store);
	esint insertElements(const ElementsRegionStore *store, const std::vector<char, initless_allocator<char> > &data);
	esint insertElements(const BoundaryRegionStore *store, const std::vector<char, initless_allocator<char> > &data);

	void insertData(const step::Step &step, NamedData *data, esint nindices, esint *indices);
	void insertDecomposition(const ElementsRegionStore *store);
protected:
	std::string _suffix;
	VTKASCIIWritter _writer;

	std::vector<char, initless_allocator<char> > _points, _esize, _ecode;
	std::vector<std::vector<char, initless_allocator<char> > > _cells;
};

}

#endif /* SRC_OUTPUT_VISUALIZATION_COLLECTED_VTKLEGACY_H_ */
