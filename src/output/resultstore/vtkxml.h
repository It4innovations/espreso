
#ifndef SRC_OUTPUT_RESULTSTORE_VTKXML_H_
#define SRC_OUTPUT_RESULTSTORE_VTKXML_H_

#include "../resultstore.h"

class vtkUnstructuredGrid;

namespace espreso {

namespace output {

class VTKXML: public ResultStore {

public:
	virtual ~VTKXML();
	virtual void finalize();

protected:
	VTKXML(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual void regionPreprocessing(const espreso::Region &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements);

	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data);
	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution);

	virtual void initWriter(const std::string &name, size_t points, size_t cells);
	virtual void addMesh(std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements);
	virtual void addData(size_t points, size_t cells, const DataArrays &data);
	virtual void addData(size_t points, size_t cells, const std::vector<Solution*> &solution);
	virtual void finalizeWriter();

	virtual void linkClusters(const std::string &root, const std::string &name, const DataArrays &data);
	virtual void linkClusters(const std::string &root, const std::string &name, const std::vector<Solution*> &solution, size_t points, size_t cells);

	virtual void linkSteps(const std::string &root, const std::string &name, const DataArrays &data);
	virtual void linkSteps(const std::string &root, const std::string &name, const std::vector<Solution*> &solution);

	virtual std::string format() const =0;
	virtual void store(std::ostream &os, const std::vector<eslocal> &data) =0;
	virtual void store(std::ostream &os, const std::vector<double> &data) =0;

	vtkUnstructuredGrid *_VTKGrid;
	std::ofstream *_os;
};

}
}


#endif /* SRC_OUTPUT_RESULTSTORE_VTKXML_H_ */
