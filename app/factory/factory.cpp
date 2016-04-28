
#include "factory.h"

using namespace espreso;

template<class TShape>
static void generateShape(const Options &options, Mesh *mesh)
{
	input::Settings settings(options, config::env::MPIrank ,config::env::MPIsize);

	switch (settings.shape) {
	case input::CUBE: {
		input::CubeSettings cube(options, config::env::MPIrank ,config::env::MPIsize);
		input::CubeGenerator<TShape>::load(*mesh, cube);
		break;
	}
	case input::SPHERE: {
		input::SphereSettings sphere(options, config::env::MPIrank ,config::env::MPIsize);
		input::SphereGenerator<TShape>::load(*mesh, sphere);
		break;
	}
	default: {
		ESINFO(ERROR) << "Unknown shape.";
	}
	}
}

static void generate(const Options &options, Mesh *mesh)
{
	input::Settings settings(options, config::env::MPIrank ,config::env::MPIsize);

	switch (settings.elementType) {
	case input::HEXA8: {
		generateShape<input::Hexahedron8>(options, mesh);
		break;
	}
	case input::HEXA20: {
		generateShape<input::Hexahedron20>(options, mesh);
		break;
	}
	case input::TETRA4: {
		generateShape<input::Tetrahedron4>(options, mesh);
		break;
	}
	case input::TETRA10: {
		generateShape<input::Tetrahedron10>(options, mesh);
		break;
	}
	case input::PRISMA6: {
		generateShape<input::Prisma6>(options, mesh);
		break;
	}
	case input::PRISMA15: {
		generateShape<input::Prisma15>(options, mesh);
		break;
	}
	case input::PYRAMID5: {
		generateShape<input::Pyramid5>(options, mesh);
		break;
	}
	case input::PYRAMID13: {
		generateShape<input::Pyramid13>(options, mesh);
		break;
	}
	default: {
		ESINFO(ERROR) << "Unknown element type.";
	}
	}
}


static Mesh* getMesh(const Options &options)
{
	Mesh *mesh = new Mesh();
	switch (config::mesh::input) {

	case config::mesh::ANSYS_MATSOL: {
		input::AnsysMatsol::load(*mesh, options, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::ANSYS_WORKBENCH: {
		input::AnsysWorkbench::load(*mesh, options, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::OPENFOAM: {
		input::OpenFOAM::load(*mesh, options, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::ESDATA: {
		input::Esdata::load(*mesh, options, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::GENERATOR: {
		generate(options, mesh);
		break;
	}
	}
	return mesh;
}

template<class TDiscretization>
static AssemblerBase* createAssembler(TDiscretization discretization)
{
	switch (config::assembler::physics) {

	case config::assembler::LinearElasticity: {
		return new LinearElasticity<TDiscretization>(discretization);
	}
	case config::assembler::Temperature: {
		return new Temperature<TDiscretization>(discretization);
	}
	case config::assembler::TransientElasticity: {
		return new TransientElasticity<TDiscretization>(discretization);
	}
	default:
		ESINFO(ERROR) << "Unknown assembler.";
		exit(EXIT_FAILURE);
	}
}

static AssemblerBase* getAssembler(Mesh *mesh, Mesh* &surface)
{
	switch (config::assembler::discretization) {

	case config::assembler::FEM: {
		FEM fem(*mesh);
		return createAssembler<FEM>(fem);
	}
	case config::assembler::BEM: {
		surface = new Mesh();
		mesh->getSurface(*surface);
		surface->computeFixPoints(config::mesh::fixPoints);
		BEM bem(*mesh, *surface);
		return createAssembler<BEM>(bem);
	}
	default:
		ESINFO(ERROR) << "Unknown discretization.";
		exit(EXIT_FAILURE);
	}
}

Factory::Factory(const Options &options)
:_assembler(NULL), _mesh(NULL), _surface(NULL)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);

	ESINFO(OVERVIEW) << "Run ESPRESO on " << config::env::MPIsize << " processes.";

	_mesh = getMesh(options);

	if (config::output::saveMesh) {
		output::VTK_Full::mesh(*_mesh, "mesh", config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
	}
	if (config::output::saveFixPoints) {
		output::VTK_Full::fixPoints(*_mesh, "meshFixPoints", config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
	}
	if (config::output::saveCorners) {
		output::VTK_Full::corners(*_mesh, "meshCorners", config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
	}
	if (config::output::saveDirichlet) {
		output::VTK_Full::dirichlet(*_mesh, "meshDirichlet", config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
	}
	if (config::output::saveAveraging) {
		output::VTK_Full::averaging(*_mesh, "meshAveraging", config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
	}

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

void Factory::solve()
{
	_assembler->init();

	for (int i = 0; i < config::assembler::timeSteps; i++) {
		_assembler->pre_solve_update();
		_assembler->solve(_solution);
		_assembler->post_solve_update();
	}

	_assembler->finalize();
}

void Factory::store(const std::string &file)
{
	switch (config::assembler::discretization){

	case config::assembler::FEM: {
		if (config::output::saveResults) {
			output::VTK_Full vtk(*_mesh, file);
			vtk.store(_solution, _mesh->DOFs(), config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
		}
		break;
	}

	case config::assembler::BEM: {
		if (config::output::saveResults) {
			output::VTK_Full vtk(*_surface, file);
			vtk.store(_solution, _surface->DOFs(), config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
		}
		break;
	}

	case config::assembler::API: {
		if (config::output::saveResults) {
			output::VTK_Full vtk(*_mesh, file);
			vtk.store(_solution, _mesh->DOFs(), config::output::subdomainShrinkRatio, config::output::clusterShrinkRatio);
		}
		break;
	}
	}
}


