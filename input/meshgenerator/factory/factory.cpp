
#include "factory.h"

using namespace esinput;

Generator* MeshFactory::create(int argc, char** argv, int rank, int size)
{
	FactorySettings settings(argc, argv);

	switch (settings.shape) {
	case CUBE: {
		switch (settings.elementType) {
		case HEXA8: {
			return new CubeGenerator<Hexahedron8>(argc, argv, rank, size);
			break;
		}
		case HEXA20: {
			return new CubeGenerator<Hexahedron20>(argc, argv, rank, size);
			break;
		}
		case TETRA4: {
			return new CubeGenerator<Tetrahedron4>(argc, argv, rank, size);
			break;
		}
		case TETRA10: {
			return new CubeGenerator<Tetrahedron10>(argc, argv, rank, size);
			break;
		}
		case PRISMA6: {
			return new CubeGenerator<Prisma6>(argc, argv, rank, size);
			break;
		}
		case PRISMA15: {
			return new CubeGenerator<Prisma15>(argc, argv, rank, size);
			break;
		}
		case PYRAMID5: {
			return new CubeGenerator<Pyramid5>(argc, argv, rank, size);
			break;
		}
		case PYRAMID13: {
			return new CubeGenerator<Pyramid13>(argc, argv, rank, size);
			break;
		}
		default: {
			std::cerr << "Unknown element type.\n";
			exit(EXIT_FAILURE);
		}
		}
		break;
	}
	case SPHERE: {
		switch (settings.elementType) {
		case HEXA8: {
			return new SphereGenerator<Hexahedron8>(argc, argv, rank, size);
			break;
		}
		case HEXA20: {
			return new SphereGenerator<Hexahedron20>(argc, argv, rank, size);
			break;
		}
		case TETRA4: {
			return new SphereGenerator<Tetrahedron4>(argc, argv, rank, size);
			break;
		}
		case TETRA10: {
			return new SphereGenerator<Tetrahedron10>(argc, argv, rank, size);
			break;
		}
		case PRISMA6: {
			return new SphereGenerator<Prisma6>(argc, argv, rank, size);
			break;
		}
		case PRISMA15: {
			return new SphereGenerator<Prisma15>(argc, argv, rank, size);
			break;
		}
		case PYRAMID5: {
			return new SphereGenerator<Pyramid5>(argc, argv, rank, size);
			break;
		}
		case PYRAMID13: {
			return new SphereGenerator<Pyramid13>(argc, argv, rank, size);
			break;
		}
		default: {
			std::cerr << "Unknown element type.\n";
			exit(EXIT_FAILURE);
		}
		break;
		}
		break;
	}
	default: {
		std::cerr << "Unknown shape.\n";
		exit(EXIT_FAILURE);
	}
	}
}



