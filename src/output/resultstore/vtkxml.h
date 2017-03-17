
#ifndef SRC_OUTPUT_RESULTSTORE_VTKXML_H_
#define SRC_OUTPUT_RESULTSTORE_VTKXML_H_

#include "../resultstore.h"

class vtkUnstructuredGrid;
class vtkXMLWriter;

namespace espreso {

namespace output {

class VTKXML: public ResultStore {

public:
	virtual ~VTKXML();

protected:
	VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);

	virtual void store(const std::string &name, const MeshInfo *regionInfo);

	virtual void initWriter(const std::string &name, size_t points, size_t cells);
	virtual void addMesh(const MeshInfo *regionInfo);
	virtual void addData(const DataArrays &data, const std::vector<Solution*> &solution);
	virtual void finalizeWriter();

	virtual void linkClusters(const std::string &root, const std::string &name, const MeshInfo *regionInfo);
	virtual void linkSteps(const std::string &name, const std::vector<std::pair<std::string, Step> > &steps);

	virtual std::string format() =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<std::vector<int> > &data) =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<std::vector<long> > &data) =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<std::vector<double> > &data) =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<int> &data) =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<long> &data) =0;
	virtual void storePointData(const std::string &name, size_t components, const std::vector<double> &data) =0;

	virtual void storeCellData(const std::string &name, size_t components, const std::vector<std::vector<int> > &data) =0;
	virtual void storeCellData(const std::string &name, size_t components, const std::vector<std::vector<long> > &data) =0;
	virtual void storeCellData(const std::string &name, size_t components, const std::vector<std::vector<double> > &data) =0;
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
