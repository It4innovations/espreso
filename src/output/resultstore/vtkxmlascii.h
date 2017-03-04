
#ifndef SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_
#define SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_

#include "vtkxml.h"

namespace espreso {

namespace output {

class VTKXMLASCII: public VTKXML {

public:
	VTKXMLASCII(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual ~VTKXMLASCII() {};

protected:
	virtual void store(std::ostream &os, const std::vector<eslocal> &data);
	virtual void store(std::ostream &os, const std::vector<double> &data);
	std::string format() const { return "ascii"; }
};

}
}


#endif /* SRC_OUTPUT_RESULTSTORE_VTKXMLASCII_H_ */
