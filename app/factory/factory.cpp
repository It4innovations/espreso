
#include "factory.h"

using namespace assembler;

static mesh::Mesh* getMesh(int argc, char **argv)
{
	mesh::Mesh *mesh = new mesh::Mesh();
	switch (esconfig::mesh::input) {

	case esconfig::mesh::ANSYS: {
		esinput::Ansys loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::OPENFOAM: {
		esinput::OpenFOAM loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::ESDATA_IN: {
		esinput::Esdata loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::GENERATOR: {
		esinput::MeshGenerator loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	}
	return mesh;
}

template<class TDiscretization>
static AssemblerBase* createAssembler(TDiscretization discretization)
{
	switch (esconfig::assembler::assembler) {

	case esconfig::assembler::LinearElasticity: {
		return new LinearElasticity<TDiscretization>(discretization);
	}
	case esconfig::assembler::Temperature: {
		return new LinearElasticity<TDiscretization>(discretization);
	}
	default:
		std::cerr << "Unknown assembler.\n";
		exit(EXIT_FAILURE);
	}
}

static AssemblerBase* getAssembler(mesh::Mesh *mesh, mesh::Mesh *surface, assembler::APIHolder *apiHolder)
{
	switch (esconfig::assembler::discretization) {

	case esconfig::assembler::FEM: {
		FEM fem(*mesh);
		return createAssembler<FEM>(fem);
	}
	case esconfig::assembler::BEM: {
		surface = new mesh::Mesh();
		mesh->getSurface(*surface);
		surface->computeFixPoints(esconfig::mesh::fixPoints);
		BEM bem(*mesh, *surface);
		return createAssembler<BEM>(bem);
	}
	case esconfig::assembler::API: {
		apiHolder = new APIHolder();
		esconfig::assembler::discretization = esconfig::assembler::FEM;
		assembler::AssemblerBase * assembler = getAssembler(mesh, surface, apiHolder);
		esconfig::assembler::discretization = esconfig::assembler::API;
		assembler->fillAPIHolder(apiHolder);
		API2 api(*apiHolder);
		delete assembler;
		return createAssembler<API2>(api);
	}
	default:
		std::cerr << "Unknown discretization.\n";
		exit(EXIT_FAILURE);
	}
}

Factory::Factory(int argc, char **argv)
:_assembler(NULL), _mesh(NULL), _surface(NULL), _apiHolder(NULL)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);

	_mesh = getMesh(argc, argv);
	_assembler = getAssembler(_mesh, _surface, _apiHolder);
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
	if (_apiHolder != NULL) {
		delete _apiHolder;
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

	store();
}

void Factory::store()
{
	switch (esconfig::assembler::discretization){

	case esconfig::assembler::FEM: {
		esoutput::VTK_Full vtk(*_mesh, "mesh");
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}

	case esconfig::assembler::BEM: {
		esoutput::VTK_Full vtk(*_surface, "surface");
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}

	case esconfig::assembler::API: {
		esoutput::VTK_Full vtk(*_mesh, "api");
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}
	}
}


