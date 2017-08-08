
#ifndef SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_
#define SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_

#include "vtkxml.h"

namespace espreso {

class VTKXMLASCII: public VTKXML {

public:
	VTKXMLASCII(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode = MeshInfo::EMPTY);
	virtual ~VTKXMLASCII();

protected:
	virtual std::string format() { return "ascii"; }

	virtual void storePointData(const std::string &name, size_t points, const std::vector<int> &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<long> &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<double> &data);

	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<int> &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<long> &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<double> &data);
};

}


#endif /* SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_ */
