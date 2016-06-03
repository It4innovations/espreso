
#include "factory.h"

using namespace espreso;

template<class TShape>
static void generateShape(const Configuration &configuration, Mesh *mesh)
{
	input::Settings settings(config::env::MPIrank ,config::env::MPIsize);
	ParametersReader::pickConfiguration(configuration, settings.parameters);

	switch (settings.shape) {
	case input::CUBE: {
		input::CubeSettings cube(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::CubeGenerator<TShape>::load(*mesh, cube);
		break;
	}
	case input::SPHERE: {
		input::SphereSettings sphere(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::SphereGenerator<TShape>::load(*mesh, sphere);
		break;
	}
	default: {
		ESINFO(ERROR) << "Unknown shape.";
	}
	}
}

static void generate(const Configuration &configuration, Mesh *mesh)
{
	input::Settings settings(config::env::MPIrank ,config::env::MPIsize);
	ParametersReader::pickConfiguration(configuration, settings.parameters);


	switch (settings.elementType) {
	case input::HEXA8: {
		generateShape<input::Hexahedron8>(configuration, mesh);
		break;
	}
	case input::HEXA20: {
		generateShape<input::Hexahedron20>(configuration, mesh);
		break;
	}
	case input::TETRA4: {
		generateShape<input::Tetrahedron4>(configuration, mesh);
		break;
	}
	case input::TETRA10: {
		generateShape<input::Tetrahedron10>(configuration, mesh);
		break;
	}
	case input::PRISMA6: {
		generateShape<input::Prisma6>(configuration, mesh);
		break;
	}
	case input::PRISMA15: {
		generateShape<input::Prisma15>(configuration, mesh);
		break;
	}
	case input::PYRAMID5: {
		generateShape<input::Pyramid5>(configuration, mesh);
		break;
	}
	case input::PYRAMID13: {
		generateShape<input::Pyramid13>(configuration, mesh);
		break;
	}
	default: {
		ESINFO(ERROR) << "Unknown element type.";
	}
	}
}


static Mesh* getMesh(const Configuration &configuration)
{
	Mesh *mesh = new Mesh();
	switch (config::mesh::INPUT) {

	case config::mesh::INPUTalternatives::MATSOL: {
		input::AnsysMatsol::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternatives::WORKBENCH: {
		input::AnsysWorkbench::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternatives::OPENFOAM: {
		input::OpenFOAM::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternatives::ESDATA: {
		input::Esdata::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternatives::GENERATOR: {
		generate(configuration, mesh);
		break;
	}
	default:
		ESINFO(GLOBAL_ERROR) << "Invalid user-supplied parameter: INPUT";
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
		surface->computeFixPoints(config::mesh::FIX_POINTS);
		BEM bem(*mesh, *surface);
		return createAssembler<BEM>(bem);
	}
	default:
		ESINFO(ERROR) << "Unknown discretization.";
		exit(EXIT_FAILURE);
	}
}

Factory::Factory(const Configuration &configuration)
:_assembler(NULL), _mesh(NULL), _surface(NULL)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);

	ESINFO(OVERVIEW) << "Run ESPRESO on " << config::env::MPIsize << " process(es).";
	ParametersReader::printParameters(config::parameters, config::info::verboseLevel);

	_mesh = getMesh(configuration);

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


