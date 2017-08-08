
#ifndef SRC_OUTPUT_RESULTSTORE_VTKXMLBINARY_H_
#define SRC_OUTPUT_RESULTSTORE_VTKXMLBINARY_H_

#include "vtkxml.h"

namespace espreso {

class VTKXMLBinary: public VTKXML {

public:
	VTKXMLBinary(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode = MeshInfo::EMPTY);
	virtual ~VTKXMLBinary();

protected:
	virtual std::string format() { return "binary"; }

	virtual void storePointData(const std::string &name, size_t points, const std::vector<int> &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<long> &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<double> &data);

	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<int> &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<long> &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<double> &data);
};

}

#endif /* SRC_OUTPUT_RESULTSTORE_VTKXMLBINARY_H_ */
