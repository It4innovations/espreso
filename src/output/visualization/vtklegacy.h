
#ifndef SRC_OUTPUT_VISUALIZATION_COLLECTED_VTKLEGACY_H_
#define SRC_OUTPUT_VISUALIZATION_COLLECTED_VTKLEGACY_H_

#include "basis/containers/allocators.h"
#include "visualization.h"
#include "writer/vtkwritter.h"
#include <string>

namespace espreso {

class RegionStore;
class ElementsRegionStore;
class BoundaryRegionStore;

class VTKLegacy: public Visualization {
public:
	VTKLegacy(const Mesh &mesh, bool withDecomposition = true);
	~VTKLegacy();

	void updateMesh();
	void updateSolution();

protected:
	void insertHeader();
	void insertPoints(const RegionStore *store);
	esint insertElements(const ElementsRegionStore *store, const std::vector<char, initless_allocator<char> > &data);
	esint insertElements(const BoundaryRegionStore *store, const std::vector<char, initless_allocator<char> > &data);

	void insertData(NamedData *data, esint nindices, esint *indices);
	void insertDecomposition(const ElementsRegionStore *store);
protected:
	bool _withDecomposition;
	std::string _suffix;
	VTKASCIIWritter _writer;

	std::vector<char, initless_allocator<char> > _points, _esize, _ecode;
	std::vector<std::vector<char, initless_allocator<char> > > _cells;
};

}

#endif /* SRC_OUTPUT_VISUALIZATION_COLLECTED_VTKLEGACY_H_ */
