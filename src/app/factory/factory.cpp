
#include "factory.h"
#include "esbasis.h"
#include "esinput.h"
#include "esoutput.h"

namespace espreso {

template<class TShape>
static void generateShape(Mesh &mesh, const Configuration &configuration, GeneratorShape shape)
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
	default: {
		ESINFO(ERROR) << "Unknown shape.";
	}
	}
}

static void generate(Mesh &mesh, const Configuration &configuration, ElementType eType, GeneratorShape shape)
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

void Factory::readParameters(const Configuration &configuration)
{
	std::vector<Parameter> parameters = {
		{ "SHAPE", shape      , "Generated shape.", {
				{ "CUBE"  , GeneratorShape::CUBE  , "Cubic." },
				{ "SPHERE", GeneratorShape::SPHERE, "Spherical." },
				{ "PLANE" , GeneratorShape::PLANE , "2D plane." }
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
				{ "LINEAR_ELASTICITY", PhysicsAssembler::LINEAR_ELASTICITY, "Linear elasticity." },
				{ "TEMPERATURE", PhysicsAssembler::TEMPERATURE, "Temperature." },
				{ "TRANSIENT_ELASTICITY", PhysicsAssembler::TRANSIENT_ELASTICITY, "Transient elasticity." },
				{ "ADVECTION_DIFFUSION", PhysicsAssembler::ADVECTION_DIFFUSION, "Advection diffusion"}
		} }
	};

	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

static void fillMesh(const Configuration &configuration, Mesh &mesh, ElementType eType, GeneratorShape shape)
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

Factory::Factory(const Configuration &configuration)
{
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

	if (config::output::SAVE_MESH) {
		std::stringstream ss;
		ss << "mesh" << config::env::MPIrank;
		output::VTK_Full::mesh(mesh, ss.str(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_FIX_POINTS) {
		std::stringstream ss;
		ss << "meshFixPoints" << config::env::MPIrank;
		output::VTK_Full::fixPoints(mesh, ss.str(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_CORNERS) {
		std::stringstream ss;
		ss << "meshCorners" << config::env::MPIrank;
		output::VTK_Full::corners(mesh, ss.str(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_DIRICHLET) {
		std::stringstream ss;
		ss << "meshDirichlet" << config::env::MPIrank;
		output::VTK_Full::dirichlet(mesh, ss.str(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}
	if (config::output::SAVE_AVERAGING) {
		std::stringstream ss;
		ss << "meshAveraging" << config::env::MPIrank;
		output::VTK_Full::averaging(mesh, ss.str(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
	}

	switch (physics) {
	case PhysicsAssembler::LINEAR_ELASTICITY:
		instance = new LinearInstance<EqualityConstraints, LinearElasticity>(mesh);
		break;
	case PhysicsAssembler::TEMPERATURE:
		instance = new LinearInstance<EqualityConstraints, Temperature>(mesh);
		break;
	case PhysicsAssembler::TRANSIENT_ELASTICITY:
		instance = new DynamicsInstance<EqualityConstraints, TransientElasticity>(mesh);
		break;
	case PhysicsAssembler::ADVECTION_DIFFUSION:
		instance = new LinearInstance<EqualityConstraints, AdvectionDiffusion>(mesh);
		break;
	}
}

void Factory::solve(const std::string &outputFile)
{
	instance->init();

	for (int i = 0; i < config::solver::TIME_STEPS; i++) {
		instance->pre_solve_update(_solution);
		instance->solve(_solution);
		instance->post_solve_update(_solution);

		if (config::output::SAVE_RESULTS) {
			std::stringstream ss;
			ss << outputFile << config::env::MPIrank;
			if (config::solver::TIME_STEPS > 1) {
				ss << "_" << i;
			}

			output::VTK_Full vtk(mesh, ss.str());
			vtk.store(_solution, config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);
		}
	}
	instance->finalize();
}

}


