
#include "stl.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/surfacestore.h"

using namespace espreso;

STL::STL(const Mesh &mesh)
: Visualization(mesh)
{

}

STL::~STL()
{

}

void STL::updateMesh()
{
	std::string filename = _directory + "file.stl";
	std::string name = _path + filename;

	int size, gsize;
	size = _mesh.surface->triangles->structures();
	Communication::reduce(&size, &gsize, 1, MPI_INT, MPI_SUM, 0, MPITools::asynchronous);

	if (info::mpi::rank == 0) {
		_writer.storeHeader("surface");
		_writer.storeSize(gsize);
	}

	for (auto t = _mesh.surface->triangles->cbegin(); t != _mesh.surface->triangles->cend(); ++t) {
		Point p[3] = {
				_mesh.nodes->coordinates->datatarray()[t->at(0)],
				_mesh.nodes->coordinates->datatarray()[t->at(1)],
				_mesh.nodes->coordinates->datatarray()[t->at(2)],
		};
		Point n = Point::cross(p[1] - p[0], p[2] - p[0]);

		_writer.beginFace(n.x, n.y, n.z);
		_writer.addVertex(p[0].x, p[0].y, p[0].z);
		_writer.addVertex(p[1].x, p[1].y, p[1].z);
		_writer.addVertex(p[2].x, p[2].y, p[2].z);
		_writer.endFace();
	}

	_writer.groupData();

	if (info::mpi::rank + 1 == info::mpi::size) {
		_writer.storeFooter("surface");
	}

	_writer.commitFile(name);
	_writer.reorder();
	_writer.write();
}


void STL::updateSolution()
{
	// TODO
}
