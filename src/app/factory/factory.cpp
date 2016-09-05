
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
	case GeneratorShape::CUBES: {
		input::CubesSettings cubes(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::CubesGenerator<TShape>::load(mesh, cubes);
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

void Factory::readParameters(const Configuration &configuration)
{
	std::vector<Parameter> parameters = {
		{ "SHAPE", shape      , "Generated shape.", {
				{ "CUBE"  , GeneratorShape::CUBE  , "Cubic." },
				{ "SPHERE", GeneratorShape::SPHERE, "Spherical." },
				{ "PLANE" , GeneratorShape::PLANE , "2D plane." },
				{ "CUBES" , GeneratorShape::CUBES , "Cubes with Mortar interface." }
		} },
		{ "ELEMENT_TYPE", eType, "The type of generated element.", {
				{ "HEXA8"    , ElementType::HEXA8    , "Hexahedron."},
				{ "HEXA20"   , ElementType::HEXA20   , "Hexahedron with midpoints."},
				{ "TETRA4"   , ElementType::TETRA4   , "Tetrahedron."},
				{ "TETRA10"  , ElementType::TETRA10  , "Tetrahedron with midpoints."},
				{ "PRISMA6"  , ElementType::PRISMA6  , "Prisma."},
				{ "PRISMA15" , ElementType::PRISMA15 , "Prisma with midpoints."},
				{ "PYRAMID5" , ElementType::PYRAMID5 , "Pyramid."},
				{ "PYRAMID13", ElementType::PYRAMID13, "Pyramid with midpoints."},

				{ "SQUARE4"  , ElementType::SQUARE4  , "Square."},
				{ "SQUARE8"  , ElementType::SQUARE8  , "Square with midpoints."},
				{ "TRIANGLE3", ElementType::TRIANGLE3, "Triangle."},
				{ "TRIANGLE6", ElementType::TRIANGLE6, "Triangle with midpoints."},
		} },
		{ "PHYSICS", physics, "Physics used for compose matrices", {
				{ "LINEAR_ELASTICITY_2D", PhysicsAssembler::LINEAR_ELASTICITY_2D, "2D linear elasticity." },
				{ "LINEAR_ELASTICITY_3D", PhysicsAssembler::LINEAR_ELASTICITY_3D, "3D linear elasticity." },
				{ "TRANSIENT_ELASTICITY_2D", PhysicsAssembler::TRANSIENT_ELASTICITY_2D, "2D transient elasticity." },
				{ "TRANSIENT_ELASTICITY_3D", PhysicsAssembler::TRANSIENT_ELASTICITY_3D, "3D transient elasticity." },
				{ "ADVECTION_DIFFUSION_2D", PhysicsAssembler::ADVECTION_DIFFUSION_2D, "2D advection diffusion"},
				{ "ADVECTION_DIFFUSION_3D", PhysicsAssembler::ADVECTION_DIFFUSION_3D, "3D advection diffusion"},
				{ "STOKES", PhysicsAssembler::STOKES, "Stokes"}
		} }
	};

	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
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

template<class TConstraints, class TPhysics, class TInstance>
static void createInstance(Instance* &instance, Mesh &mesh)
{
	switch (config::solver::FETI_METHOD) {
	case config::solver::FETI_METHODalternative::TOTAL_FETI:
	case config::solver::FETI_METHODalternative::HYBRID_FETI:
		instance = new TInstance(mesh);
		break;
	case config::solver::FETI_METHODalternative::HYPRE:
		instance = new HypreInstance<TConstraints, TPhysics>(mesh);
		break;
	}
}

Factory::Factory(const Configuration &configuration)
{
	physics = PhysicsAssembler::LINEAR_ELASTICITY_3D;

	readParameters(configuration);

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

	if (config::output::SAVE_MESH) {
		output::VTK::mesh(mesh, "mesh", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_FIX_POINTS) {
		output::VTK::fixPoints(mesh, "meshFixPoints", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_CORNERS) {
		output::VTK::corners(mesh, "meshCorners", config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_PROPERTIES) {
//		output::VTK::properties(mesh, "meshDirichletX",
//				{ Property::DISPLACEMENT_X },
//				config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
//
//		output::VTK::properties(mesh, "meshDirichletY",
//				{ Property::DISPLACEMENT_Y },
//				config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
//
//		output::VTK::properties(mesh, "meshDirichletZ",
//				{ Property::DISPLACEMENT_Z },
//				config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);

		output::VTK::properties(mesh, "meshTranslationMotion",
				{ Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y },
				config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}

	switch (physics) {
	case PhysicsAssembler::LINEAR_ELASTICITY_2D:
		createInstance<EqualityConstraints, LinearElasticity2D, LinearInstance<EqualityConstraints, LinearElasticity2D> >(instance, mesh);
		break;
	case PhysicsAssembler::LINEAR_ELASTICITY_3D:
		createInstance<EqualityConstraints, LinearElasticity3D, LinearInstance<EqualityConstraints, LinearElasticity3D> >(instance, mesh);
		break;
	case PhysicsAssembler::TRANSIENT_ELASTICITY_3D:
		createInstance<EqualityConstraints, TransientElasticity, LinearInstance<EqualityConstraints, TransientElasticity> >(instance, mesh);
		break;
	case PhysicsAssembler::ADVECTION_DIFFUSION_2D:
		createInstance<EqualityConstraints, AdvectionDiffusion2D, LinearInstance<EqualityConstraints, AdvectionDiffusion2D> >(instance, mesh);
		break;
	case PhysicsAssembler::STOKES:
		createInstance<EqualityConstraints, Stokes, LinearInstance<EqualityConstraints, Stokes> >(instance, mesh);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown Physics";
	}
}

void Factory::solve()
{
	output::Store *results = NULL;

	instance->init();
	if (config::output::SAVE_RESULTS) {
		results = new output::VTK(mesh, "results");
	}

	for (int i = 0; i < config::solver::TIME_STEPS; i++) {
		instance->pre_solve_update(_solution);
		instance->solve(_solution);
		instance->post_solve_update(_solution);
		if (config::output::SAVE_RESULTS) {
			results->store(_solution, config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
		}
	}

	instance->finalize();
	if (results != NULL) {
		delete results;
	}
}

double Factory::norm() const
{
	double n = 0, sum = 0;
	for (size_t i = 0; i < _solution.size(); i++) {
		for (size_t j = 0; j < _solution[i].size(); j++) {
			n += _solution[i][j] * _solution[i][j];
		}
	}

	MPI_Allreduce(&n, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(sum);
}

}


