
#include "factory.h"

using namespace espreso;

template<class TShape>
static void generateShape(const Configuration &configuration, Mesh *mesh)
{
	input::Settings settings(config::env::MPIrank ,config::env::MPIsize);
	ParametersReader::pickConfiguration(configuration, settings.parameters);

	switch (settings.shape) {
	case input::GeneratorShape::CUBE: {
		input::CubeSettings cube(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::CubeGenerator<TShape>::load(*mesh, cube);
		break;
	}
	case input::GeneratorShape::SPHERE: {
		input::SphereSettings sphere(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::SphereGenerator<TShape>::load(*mesh, sphere);
		break;
	}
	case input::GeneratorShape::PLANE: {
		input::PlaneSettings plane(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::PlaneGenerator<TShape>::load(*mesh, plane);
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


	switch (settings.eType) {
	case input::ElementType::HEXA8: {
		generateShape<input::Hexahedron8>(configuration, mesh);
		break;
	}
	case input::ElementType::HEXA20: {
		generateShape<input::Hexahedron20>(configuration, mesh);
		break;
	}
	case input::ElementType::TETRA4: {
		generateShape<input::Tetrahedron4>(configuration, mesh);
		break;
	}
	case input::ElementType::TETRA10: {
		generateShape<input::Tetrahedron10>(configuration, mesh);
		break;
	}
	case input::ElementType::PRISMA6: {
		generateShape<input::Prisma6>(configuration, mesh);
		break;
	}
	case input::ElementType::PRISMA15: {
		generateShape<input::Prisma15>(configuration, mesh);
		break;
	}
	case input::ElementType::PYRAMID5: {
		generateShape<input::Pyramid5>(configuration, mesh);
		break;
	}
	case input::ElementType::PYRAMID13: {
		generateShape<input::Pyramid13>(configuration, mesh);
		break;
	}
	case input::ElementType::SQUARE4: {
		generateShape<input::Square4>(configuration, mesh);
		break;
	}
	case input::ElementType::SQUARE8: {
		generateShape<input::Square8>(configuration, mesh);
		break;
	}
	case input::ElementType::TRIANGLE3: {
		generateShape<input::Triangle3>(configuration, mesh);
		break;
	}
	case input::ElementType::TRIANGLE6: {
		generateShape<input::Triangle6>(configuration, mesh);
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

	case config::mesh::INPUTalternative::MATSOL: {
		input::AnsysMatsol::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternative::WORKBENCH: {
		input::AnsysWorkbench::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternative::OPENFOAM: {
		input::OpenFOAM::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternative::ESDATA: {
		input::Esdata::load(*mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	}
	case config::mesh::INPUTalternative::GENERATOR: {
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
	switch (config::assembler::PHYSICS) {

	case config::assembler::PHYSICSalternative::LINEAR_ELASTICITY: {
		return new LinearElasticity<TDiscretization>(discretization);
	}
	case config::assembler::PHYSICSalternative::TEMPERATURE: {
		return new Temperature<TDiscretization>(discretization);
	}
	case config::assembler::PHYSICSalternative::TRANSIENT_ELASTICITY: {
		return new TransientElasticity<TDiscretization>(discretization);
	}
	default:
		ESINFO(ERROR) << "Unknown assembler.";
		exit(EXIT_FAILURE);
	}
}

static AssemblerBase* getAssembler(Mesh *mesh, Mesh* &surface)
{
	switch (config::assembler::DISCRETIZATION) {

	case config::assembler::DISCRETIZATIONalternative::FEM: {
		FEM fem(*mesh);
		return createAssembler<FEM>(fem);
	}
	case config::assembler::DISCRETIZATIONalternative::BEM: {
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

	_mesh = getMesh(configuration);

	if (config::output::SAVE_MESH) {
		output::VTK_Full::mesh(*_mesh, "mesh", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_FIX_POINTS) {
		output::VTK_Full::fixPoints(*_mesh, "meshFixPoints", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_CORNERS) {
		output::VTK_Full::corners(*_mesh, "meshCorners", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_DIRICHLET) {
		output::VTK_Full::dirichlet(*_mesh, "meshDirichlet", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_AVERAGING) {
		output::VTK_Full::averaging(*_mesh, "meshAveraging", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
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

	for (int i = 0; i < config::assembler::TIME_STEPS; i++) {
		_assembler->pre_solve_update();
		_assembler->solve(_solution);
		_assembler->post_solve_update();
	}

	_assembler->finalize();
}

void Factory::store(const std::string &file)
{
	switch (config::assembler::DISCRETIZATION){

	case config::assembler::DISCRETIZATIONalternative::FEM: {
		if (config::output::SAVE_RESULTS) {
			output::VTK_Full vtk(*_mesh, file);
			vtk.store(_solution, _mesh->DOFs(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
		}
		break;
	}

	case config::assembler::DISCRETIZATIONalternative::BEM: {
		if (config::output::SAVE_RESULTS) {
			output::VTK_Full vtk(*_surface, file);
			vtk.store(_solution, _surface->DOFs(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
		}
		break;
	}

	case config::assembler::DISCRETIZATIONalternative::API: {
		if (config::output::SAVE_RESULTS) {
			output::VTK_Full vtk(*_mesh, file);
			vtk.store(_solution, _mesh->DOFs(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
		}
		break;
	}
	}
}


