
#ifndef SRC_OUTPUT_RESULTSTORE_VTKXML_H_
#define SRC_OUTPUT_RESULTSTORE_VTKXML_H_

#include "../resultstore.h"

class vtkUnstructuredGrid;
class vtkXMLWriter;

namespace espreso {
namespace output {

struct DataArrays;

class VTKXML: public ResultStore {

public:
	virtual ~VTKXML();

protected:
	VTKXML(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode = MeshInfo::EMPTY);

	virtual std::string store(const std::string &name, const RegionData &regionData);

	virtual std::string initWriter(const std::string &name, size_t points, size_t cells);
	virtual void addMesh(const RegionData &regionData);
	virtual void addData(const RegionData &regionData);
	virtual void finalizeWriter();

	virtual std::string linkClusters(const std::string &root, const std::string &name, const RegionData &regionData);
	virtual void linkSteps(const std::string &name, const std::vector<std::pair<Step, std::vector<std::string> > > &steps);

	virtual std::string format() =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<int> &data) =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<long> &data) =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<double> &data) =0;

	virtual void storeCellData(const std::string &name, size_t components, const std::vector<int> &data) =0;
	virtual void storeCellData(const std::string &name, size_t components, const std::vector<long> &data) =0;
	virtual void storeCellData(const std::string &name, size_t components, const std::vector<double> &data) =0;

	vtkUnstructuredGrid *_VTKGrid;
	vtkXMLWriter *_writer;
	std::vector<void*> _VTKDataArrays;

	std::ofstream *_os;
};

}
}


#endif /* SRC_OUTPUT_RESULTSTORE_VTKXML_H_ */
