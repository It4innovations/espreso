
#include "factory.h"

using namespace assembler;

template<class TShape>
static void generateShape(const Options &options, mesh::Mesh *mesh)
{
	esinput::Settings settings(options, esconfig::MPIrank ,esconfig::MPIsize);

	switch (settings.shape) {
	case esinput::CUBE: {
		esinput::CubeSettings cube(options, esconfig::MPIrank ,esconfig::MPIsize);
		esinput::CubeGenerator<TShape> generator(cube);
		generator.load(*mesh);
		break;
	}
	case esinput::SPHERE: {
		esinput::SphereSettings sphere(options, esconfig::MPIrank ,esconfig::MPIsize);
		esinput::SphereGenerator<TShape> generator(sphere);
		generator.load(*mesh);
		break;
	}
	default: {
		ESLOG(eslog::ERROR) << "Unknown shape.";
	}
	}
}

static void generate(const Options &options, mesh::Mesh *mesh)
{
	esinput::Settings settings(options, esconfig::MPIrank ,esconfig::MPIsize);

	switch (settings.elementType) {
	case esinput::HEXA8: {
		generateShape<esinput::Hexahedron8>(options, mesh);
		break;
	}
	case esinput::HEXA20: {
		generateShape<esinput::Hexahedron20>(options, mesh);
		break;
	}
	case esinput::TETRA4: {
		generateShape<esinput::Tetrahedron4>(options, mesh);
		break;
	}
	case esinput::TETRA10: {
		generateShape<esinput::Tetrahedron10>(options, mesh);
		break;
	}
	case esinput::PRISMA6: {
		generateShape<esinput::Prisma6>(options, mesh);
		break;
	}
	case esinput::PRISMA15: {
		generateShape<esinput::Prisma15>(options, mesh);
		break;
	}
	case esinput::PYRAMID5: {
		generateShape<esinput::Pyramid5>(options, mesh);
		break;
	}
	case esinput::PYRAMID13: {
		generateShape<esinput::Pyramid13>(options, mesh);
		break;
	}
	default: {
		ESLOG(eslog::ERROR) << "Unknown element type.";
	}
	}
}



static mesh::Mesh* getMesh(const Options &options)
{
	mesh::Mesh *mesh = new mesh::Mesh();
	switch (esconfig::mesh::input) {

	case esconfig::mesh::ANSYS_MATSOL: {
		esinput::AnsysMatsol loader(options, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::ANSYS_WORKBENCH: {
		esinput::AnsysWorkbench loader(options, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::OPENFOAM: {
		esinput::OpenFOAM loader(options, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::ESDATA_IN: {
		esinput::Esdata loader(options, esconfig::MPIrank, esconfig::MPIsize);
		loader.load(*mesh);
		break;
	}
	case esconfig::mesh::GENERATOR: {
		generate(options, mesh);
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
		ESLOG(eslog::ERROR) << "Unknown assembler.";
		exit(EXIT_FAILURE);
	}
}

static AssemblerBase* getAssembler(mesh::Mesh *mesh, mesh::Mesh *surface)
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
	default:
		ESLOG(eslog::ERROR) << "Unknown discretization.";
		exit(EXIT_FAILURE);
	}
}

Factory::Factory(const Options &options)
:_assembler(NULL), _mesh(NULL), _surface(NULL)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);

	_mesh = getMesh(options);
	_assembler = getAssembler(_mesh, _surface);
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


