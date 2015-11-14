
#include "factory.h"

static void load(mesh::Mesh *mesh, int argc, char **argv)
{
	switch (esconfig::mesh::input) {

	case esconfig::mesh::ANSYS: {
		esinput::Ansys loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		mesh->partitiate(esconfig::mesh::subdomains, esconfig::mesh::fixPoints);
		break;
	}
	case esconfig::mesh::OPENFOAM: {
		esinput::OpenFOAM loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		mesh->partitiate(esconfig::mesh::subdomains, esconfig::mesh::fixPoints);
		break;
	}
	case esconfig::mesh::ESDATA_IN: {
		esinput::Esdata loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		mesh->partitiate(esconfig::mesh::subdomains, esconfig::mesh::fixPoints);
		break;
	}
	case esconfig::mesh::GENERATOR: {
		esinput::MeshGenerator loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}

	}
}

using namespace assembler;

Factory::Factory(int argc, char **argv)
:_assembler(NULL), _mesh(NULL), _surface(NULL)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);

	_mesh = new mesh::Mesh(esconfig::MPIrank, esconfig::MPIsize);
	load(_mesh, argc, argv);

	switch (esconfig::assembler::assembler) {

	case esconfig::assembler::LinearElasticity: {
		switch (esconfig::assembler::discretization) {

		case esconfig::assembler::FEM: {
			FEM fem(*_mesh);
			_assembler = new LinearElasticity<FEM>(fem);
			break;
		}
		case esconfig::assembler::BEM: {
			_surface = new mesh::SurfaceMesh(*_mesh);
			_surface->computeFixPoints(esconfig::mesh::fixPoints);
			BEM bem(*_mesh, *_surface);
			_assembler = new LinearElasticity<BEM>(bem);
			break;
		}
		}
		break;
	}

	case esconfig::assembler::Temperature: {
		switch (esconfig::assembler::discretization) {

		case esconfig::assembler::FEM: {
			FEM fem(*_mesh);
			_assembler = new Temperature<FEM>(fem);
			break;
		}
		case esconfig::assembler::BEM: {
			_surface = new mesh::SurfaceMesh(*_mesh);
			_surface->computeFixPoints(esconfig::mesh::fixPoints);
			BEM bem(*_mesh, *_surface);
			_assembler = new Temperature<BEM>(bem);
			break;
		}
		}
		break;
	}
	}

}

Factory::~Factory()
{
	delete _assembler;
	if (_mesh != NULL) {
		delete _mesh;
	}
	if (_surface != NULL) {
		delete _surface;
	}
}

void Factory::solve(eslocal steps)
{
	_assembler->init();

	for (int i = 0; i < steps; i++) {
		_assembler->pre_solve_update();
		_assembler->solve(_solution);
		_assembler->post_solve_update();
	}

	_assembler->finalize();
}

void Factory::store(const char *file)
{
	switch (esconfig::assembler::discretization){

	case esconfig::assembler::FEM: {
		esoutput::VTK_Full vtk(*_mesh, file);
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}

	case esconfig::assembler::BEM: {
		esoutput::VTK_Full vtk(*_surface, "surface");
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}
	}
}


