
#ifndef SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_
#define SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_

#include "vtkxml.h"

namespace espreso {

namespace output {

class VTKXMLASCII: public VTKXML {

public:
	VTKXMLASCII(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual ~VTKXMLASCII();

protected:
	virtual std::string format() { return "ascii"; }

	virtual void storePointData(const std::string &name, size_t points, const std::vector<std::vector<int> > &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<std::vector<long> > &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<std::vector<double> > &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<int> &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<long> &data);
	virtual void storePointData(const std::string &name, size_t points, const std::vector<double> &data);

	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<int> > &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<long> > &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<std::vector<double> > &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<int> &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<long> &data);
	virtual void storeCellData(const std::string &name, size_t cells, const std::vector<double> &data);
};

}
}


#endif /* SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_ */
