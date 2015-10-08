
#ifndef OUTPUT_VTK_VTK_H_
#define OUTPUT_VTK_VTK_H_

#include "../store.h"

namespace esoutput {

class VTK: public ResultStore {

public:
	VTK(const std::string &path, int rank, int size): _file(path), _rank(rank), _size(size) { };

	void store(const mesh::Mesh &mesh, double shrinkSubdomain, double shringCluster);
	void store(const mesh::Mesh &mesh, std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shringCluster);

private:
	void head();
	void coordinates(const mesh::Coordinates &coordinates, double shrinkSubdomain, double shringCluster);
	void elements(const mesh::Mesh &mesh);
	void coordinatesDisplacement(const std::vector<std::vector<double> > &displacement);

	std::ofstream _vtk;
	std::string _file;
	int _rank;
	int _size;
};


}


#endif /* OUTPUT_VTK_VTK_H_ */
