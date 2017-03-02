
#ifndef SRC_OUTPUT2_RESULTSTORE_VTKXMLBINARY_H_
#define SRC_OUTPUT2_RESULTSTORE_VTKXMLBINARY_H_

#include "vtkxml.h"

namespace espreso {

namespace output {

class VTKXMLBinary: public VTKXML {

public:
	VTKXMLBinary(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual ~VTKXMLBinary() {};

protected:
	virtual void store(std::ostream &os, const std::vector<eslocal> &data);
	virtual void store(std::ostream &os, const std::vector<double> &data);
	std::string format() const { return "binary"; }
};

}
}

#endif /* SRC_OUTPUT2_RESULTSTORE_VTKXMLBINARY_H_ */
