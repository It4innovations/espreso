
#include "factory.h"
#include "esbasis.h"
#include "esinput.h"

namespace espreso {

template<class TShape>
static void generateShape(Mesh &mesh, const ArgsConfiguration &configuration, GeneratorShape shape)
{
	switch (shape) {
	case GeneratorShape::CUBE: {
		input::CubeSettings cube(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::CubeGenerator<TShape>::load(mesh, cube);
		break;
	}
	case GeneratorShape::SPHERE: {
		input::SphereSettings sphere(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::SphereGenerator<TShape>::load(mesh, sphere);
		break;
	}
	case GeneratorShape::PLANE: {
		input::PlaneSettings plane(configuration, config::env::MPIrank ,config::env::MPIsize);
		input::PlaneGenerator<TShape>::load(mesh, plane);
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

static void generate(Mesh &mesh, const ArgsConfiguration &configuration, ElementType eType, GeneratorShape shape)
{
	switch (eType) {
	case ElementType::HEXA8:
		generateShape<input::Hexahedron8>(mesh, configuration, shape);
		break;
	case ElementType::HEXA20:
		generateShape<input::Hexahedron20>(mesh, configuration, shape);
		break;
	case ElementType::TETRA4:
		generateShape<input::Tetrahedron4>(mesh, configuration, shape);
		break;
	case ElementType::TETRA10:
		generateShape<input::Tetrahedron10>(mesh, configuration, shape);
		break;
	case ElementType::PRISMA6:
		generateShape<input::Prisma6>(mesh, configuration, shape);
		break;
	case ElementType::PRISMA15:
		generateShape<input::Prisma15>(mesh, configuration, shape);
		break;
	case ElementType::PYRAMID5:
		generateShape<input::Pyramid5>(mesh, configuration, shape);
		break;
	case ElementType::PYRAMID13:
		generateShape<input::Pyramid13>(mesh, configuration, shape);
		break;
	case ElementType::SQUARE4:
		generateShape<input::Square4>(mesh, configuration, shape);
		break;
	case ElementType::SQUARE8:
		generateShape<input::Square8>(mesh, configuration, shape);
		break;
	case ElementType::TRIANGLE3:
		generateShape<input::Triangle3>(mesh, configuration, shape);
		break;
	case ElementType::TRIANGLE6:
		generateShape<input::Triangle6>(mesh, configuration, shape);
		break;
	default:
		ESINFO(ERROR) << "Unknown element type.";
	}
}

void Factory::readParameters(const ArgsConfiguration &configuration)
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

	if (config::mesh::INPUT == config::mesh::INPUTalternative::GENERATOR) {
		ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
	}
}

static void fillMesh(const ArgsConfiguration &configuration, Mesh &mesh, ElementType eType, GeneratorShape shape)
{
	switch (config::mesh::INPUT) {
	case config::mesh::INPUTalternative::MATSOL:
		input::AnsysMatsol::load(mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	case config::mesh::INPUTalternative::WORKBENCH:
		input::AnsysWorkbench::load(mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	case config::mesh::INPUTalternative::OPENFOAM:
		input::OpenFOAM::load(mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	case config::mesh::INPUTalternative::ESDATA:
		input::Esdata::load(mesh, configuration, config::env::MPIrank, config::env::MPIsize);
		break;
	case config::mesh::INPUTalternative::GENERATOR:
		generate(mesh, configuration, eType, shape);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Invalid user-supplied parameter: INPUT";
	}
}

template<class TPhysics, class TInstance>
static void createInstance(Instance* &instance, Mesh &mesh)
{
	switch (config::solver::FETI_METHOD) {
	case config::solver::FETI_METHODalternative::TOTAL_FETI:
	case config::solver::FETI_METHODalternative::HYBRID_FETI:
		instance = new TInstance(mesh);
		break;
	case config::solver::FETI_METHODalternative::HYPRE:
		instance = new HypreInstance<TPhysics>(mesh);
		break;
	}
}

Factory::Factory(const ArgsConfiguration &configuration)
{
	physics = PhysicsAssembler::LINEAR_ELASTICITY_3D;

	readParameters(configuration);

	switch (config::assembler::DISCRETIZATION) {
	case config::assembler::DISCRETIZATIONalternative::FEM:
		fillMesh(configuration, mesh, eType, shape);
		break;
	case config::assembler::DISCRETIZATIONalternative::BEM: {
		Mesh fullMesh;
		fillMesh(configuration, fullMesh, eType, shape);
		fullMesh.getSurface(mesh);
		} break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown discretization";
	}

	switch (physics) {
	case PhysicsAssembler::LINEAR_ELASTICITY_2D:
		createInstance<LinearElasticity2D, LinearInstance<LinearElasticity2D> >(instance, mesh);
		break;
	case PhysicsAssembler::LINEAR_ELASTICITY_3D:
		createInstance<LinearElasticity3D, LinearInstance<LinearElasticity3D> >(instance, mesh);
		break;
	case PhysicsAssembler::TRANSIENT_ELASTICITY_3D:
		createInstance<TransientElasticity, DynamicsInstance<TransientElasticity> >(instance, mesh);
		break;
	case PhysicsAssembler::ADVECTION_DIFFUSION_2D:
		createInstance<AdvectionDiffusion2D, LinearInstance<AdvectionDiffusion2D> >(instance, mesh);
		break;
	case PhysicsAssembler::ADVECTION_DIFFUSION_3D:
		createInstance<AdvectionDiffusion3D, LinearInstance<AdvectionDiffusion3D> >(instance, mesh);
		break;
	case PhysicsAssembler::STOKES:
		createInstance<Stokes, LinearInstance<Stokes> >(instance, mesh);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown Physics";
	}
}

void Factory::solve(const std::string &outputFile)
{
	instance->init();

	for (int i = 0; i < config::solver::TIME_STEPS; i++) {
		instance->solve(_solution);
	}

	instance->finalize();
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


