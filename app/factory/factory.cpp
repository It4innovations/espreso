
#include "factory.h"

using namespace assembler;

template<class TShape>
static void generateShape(int argc, char **argv, mesh::Mesh *mesh)
{
	esinput::Settings settings(argc, argv, 0 ,1);

	switch (settings.shape) {
	case esinput::CUBE: {
		esinput::CubeSettings cube(argc, argv, 0, 1);
		esinput::CubeGenerator<TShape> generator(cube);
		generator.load(*mesh);
		break;
	}
	case esinput::SPHERE: {
		esinput::SphereSettings sphere(argc, argv, 0 ,1);
		esinput::SphereGenerator<TShape> generator(sphere);
		generator.load(*mesh);
		break;
	}
	default: {
		std::cerr << "Unknown shape.\n";
		exit(EXIT_FAILURE);
	}
	}
}

static void generate(int argc, char **argv, mesh::Mesh *mesh)
{
	esinput::Settings settings(argc, argv, 0 ,1);

	switch (settings.elementType) {
	case esinput::HEXA8: {
		generateShape<esinput::Hexahedron8>(argc, argv, mesh);
		break;
	}
	case esinput::HEXA20: {
		generateShape<esinput::Hexahedron20>(argc, argv, mesh);
		break;
	}
	case esinput::TETRA4: {
		generateShape<esinput::Tetrahedron4>(argc, argv, mesh);
		break;
	}
	case esinput::TETRA10: {
		generateShape<esinput::Tetrahedron10>(argc, argv, mesh);
		break;
	}
	case esinput::PRISMA6: {
		generateShape<esinput::Prisma6>(argc, argv, mesh);
		break;
	}
	case esinput::PRISMA15: {
		generateShape<esinput::Prisma15>(argc, argv, mesh);
		break;
	}
	case esinput::PYRAMID5: {
		generateShape<esinput::Pyramid5>(argc, argv, mesh);
		break;
	}
	case esinput::PYRAMID13: {
		generateShape<esinput::Pyramid13>(argc, argv, mesh);
		break;
	}
	default: {
		std::cerr << "Unknown element type.\n";
		exit(EXIT_FAILURE);
	}
	}
}



static mesh::Mesh* getMesh(int argc, char **argv)
{
	mesh::Mesh *mesh = new mesh::Mesh();
	switch (esconfig::mesh::input) {

	case esconfig::mesh::ANSYS_MATSOL: {
		esinput::AnsysMatsol loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::ANSYS_WORKBENCH: {
		esinput::AnsysWorkbench loader(argc, argv, esconfig::MPIrank, esconfig::MPIsize);
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
		generate(argc, argv, mesh);
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
		API api(*apiHolder);
		delete assembler;
		return createAssembler<API>(api);
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
}

void Factory::store(const std::string &file)
{
	switch (esconfig::assembler::discretization){

	case esconfig::assembler::FEM: {
		esoutput::VTK_Full vtk(*_mesh, file);
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}

	case esconfig::assembler::BEM: {
		esoutput::VTK_Full vtk(*_surface, file);
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}

	case esconfig::assembler::API: {
		esoutput::VTK_Full vtk(*_mesh, file);
		vtk.store(_solution, _assembler->DOFs(), 0.95, 0.9);
		break;
	}
	}
}


