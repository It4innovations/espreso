#include "permoncube.h"


	//	###################################################
	//	#                                                 #
	//	#             A z-coord.                          #
	//	#             |                                   #
	//	#             |            E3                     #
	//	#             |_ _ _ _ _ _ _                      #
	//	#            /     E5      /|                     #
	//	#           /_ _ _ _ _ _  / |                     #
	//	#          |      |      |  |                     #
	//	#        E4|      |      |E2|                     #
	//	#          |_ _ _ |_ _ _ |  |       y-coord.      #
	//	#          |    E1|      |  |------->             #
	//	#          |      |      | /                      #
	//	#          |_ _ _ |_ _ _ |/                       #
	//	#         /                                       #
	//	#        /       E0                               #
	//	#       /                                         #
	//	#      v  x-coord.                                #
	//	#                                                 #
	//	###################################################

void Permoncube::hexahedrons8(
		Mesh &mesh,
		Coordinates &coordinates,
		int *subdomains,
		int *elementsInSub)
{
	int nnx = subdomains[0] * elementsInSub[0] + 1;
	int nny = subdomains[1] * elementsInSub[1] + 1;
	int nnz = subdomains[2] * elementsInSub[2] + 1;
	double lenght[3] = { 30, 30, 30 };

	double stepx = lenght[0] / (nnx - 1);
	double stepy = lenght[1] / (nny - 1);
	double stepz = lenght[2] / (nnz - 1);

	idx_t index = 0;
	coordinates.resize(nnz * nny * nnx);
	for (int z = 0; z < nnz; z++) {
		for (int y = 0; y < nny; y++) {
			for (int x = 0; x < nnx; x++) {
				coordinates[index++] = Point(x * stepx, y * stepy, z * stepz);
			}
		}
	}

	idx_t indices[8];
	int offset[3];
	mesh.reserve(
			subdomains[2] * subdomains[1] * subdomains[0] *
			elementsInSub[2] * elementsInSub[1] * elementsInSub[0]
	);

	for (int subz = 0; subz < subdomains[2]; subz++) {
		for (int suby = 0; suby < subdomains[1]; suby++) {
			for (int subx = 0; subx < subdomains[0]; subx++) {
				offset[2] = subz * elementsInSub[2];
				offset[1] = suby * elementsInSub[1];
				offset[0] = subx * elementsInSub[0];
				for (int z = offset[2]; z < offset[2] + elementsInSub[2]; z++) {
					for (int y = offset[1]; y < offset[1] + elementsInSub[1]; y++) {
						for (int x = offset[0]; x < offset[0] + elementsInSub[0]; x++) {
							indices[0] = nnx * nny *  z      + nnx *  y      + x;
							indices[1] = nnx * nny *  z      + nnx *  y      + x + 1;
							indices[2] = nnx * nny *  z      + nnx * (y + 1) + x + 1;
							indices[3] = nnx * nny *  z      + nnx * (y + 1) + x;
							indices[4] = nnx * nny * (z + 1) + nnx *  y      + x;
							indices[5] = nnx * nny * (z + 1) + nnx *  y      + x + 1;
							indices[6] = nnx * nny * (z + 1) + nnx * (y + 1) + x + 1;
							indices[7] = nnx * nny * (z + 1) + nnx * (y + 1) + x;
							mesh.pushElement(new Hexahedron8(indices));
						}
					}
				}
				mesh.endPartition();
			}
		}
	}

}

void Permoncube::tetrahedrons4(
		Mesh &mesh,
		Coordinates &coordinates,
		int *subdomains,
		int *elementsInSub)
{
	int nnx = subdomains[0] * elementsInSub[0] + 1;
	int nny = subdomains[1] * elementsInSub[1] + 1;
	int nnz = subdomains[2] * elementsInSub[2] + 1;
	double lenght[3] = { 30, 30, 30 };

	double stepx = lenght[0] / (nnx - 1);
	double stepy = lenght[1] / (nny - 1);
	double stepz = lenght[2] / (nnz - 1);

	idx_t index = 0;
	coordinates.resize(nnz * nny * nnx);
	for (int z = 0; z < nnz; z++) {
		for (int y = 0; y < nny; y++) {
			for (int x = 0; x < nnx; x++) {
				coordinates[index++] = Point(x * stepx, y * stepy, z * stepz);
			}
		}
	}

	idx_t indices[8];
	idx_t tetra[5];
	int offset[3];
	mesh.reserve(
			subdomains[2] * subdomains[1] * subdomains[0] *
			elementsInSub[2] * elementsInSub[1] * elementsInSub[0] *
			6
	);

	for (int subz = 0; subz < subdomains[2]; subz++) {
		for (int suby = 0; suby < subdomains[1]; suby++) {
			for (int subx = 0; subx < subdomains[0]; subx++) {
				offset[2] = subz * elementsInSub[2];
				offset[1] = suby * elementsInSub[1];
				offset[0] = subx * elementsInSub[0];
				for (int z = offset[2]; z < offset[2] + elementsInSub[2]; z++) {
					for (int y = offset[1]; y < offset[1] + elementsInSub[1]; y++) {
						for (int x = offset[0]; x < offset[0] + elementsInSub[0]; x++) {
							indices[0] = nnx * nny *  z      + nnx *  y      + x;
							indices[1] = nnx * nny *  z      + nnx *  y      + x + 1;
							indices[2] = nnx * nny *  z      + nnx * (y + 1) + x + 1;
							indices[3] = nnx * nny *  z      + nnx * (y + 1) + x;
							indices[4] = nnx * nny * (z + 1) + nnx *  y      + x;
							indices[5] = nnx * nny * (z + 1) + nnx *  y      + x + 1;
							indices[6] = nnx * nny * (z + 1) + nnx * (y + 1) + x + 1;
							indices[7] = nnx * nny * (z + 1) + nnx * (y + 1) + x;

							tetra[0] = indices[0];
							tetra[1] = indices[2];
							tetra[2] = indices[3];
							tetra[4] = indices[4];
							mesh.pushElement(new Tetrahedron4(tetra));

							tetra[0] = indices[2];
							tetra[1] = indices[3];
							tetra[2] = indices[4];
							tetra[4] = indices[7];
							mesh.pushElement(new Tetrahedron4(tetra));

							tetra[0] = indices[6];
							tetra[1] = indices[2];
							tetra[2] = indices[4];
							tetra[4] = indices[7];
							mesh.pushElement(new Tetrahedron4(tetra));

							tetra[0] = indices[2];
							tetra[1] = indices[5];
							tetra[2] = indices[6];
							tetra[4] = indices[4];
							mesh.pushElement(new Tetrahedron4(tetra));

							tetra[0] = indices[1];
							tetra[1] = indices[5];
							tetra[2] = indices[2];
							tetra[4] = indices[4];
							mesh.pushElement(new Tetrahedron4(tetra));

							tetra[0] = indices[0];
							tetra[1] = indices[4];
							tetra[2] = indices[1];
							tetra[4] = indices[2];
							mesh.pushElement(new Tetrahedron4(tetra));
						}
					}
				}
				mesh.endPartition();
			}
		}
	}
}

void Permoncube::tetrahedrons10(
		Mesh &mesh,
		Coordinates &coordinates,
		int *subdomains,
		int *elementsInSub)
{
	int nnx = 2 * (subdomains[0] * elementsInSub[0]) + 1;
	int nny = 2 * (subdomains[1] * elementsInSub[1]) + 1;
	int nnz = 2 * (subdomains[2] * elementsInSub[2]) + 1;
	double lenght[3] = { 30, 30, 30 };

	double stepx = lenght[0] / (nnx - 1);
	double stepy = lenght[1] / (nny - 1);
	double stepz = lenght[2] / (nnz - 1);

	idx_t index = 0;
	coordinates.resize(nnz * nny * nnx);
	for (int z = 0; z < nnz; z++) {
		for (int y = 0; y < nny; y++) {
			for (int x = 0; x < nnx; x++) {
				coordinates[index++] = Point(x * stepx, y * stepy, z * stepz);
			}
		}
	}

	idx_t indices[27];
	idx_t tetra[20];
	int offset[3];
	mesh.reserve(
			subdomains[2] * subdomains[1] * subdomains[0] *
			elementsInSub[2] * elementsInSub[1] * elementsInSub[0] *
			6
	);

	for (int subz = 0; subz < subdomains[2]; subz++) {
		for (int suby = 0; suby < subdomains[1]; suby++) {
			for (int subx = 0; subx < subdomains[0]; subx++) {
				offset[2] = subz * elementsInSub[2];
				offset[1] = suby * elementsInSub[1];
				offset[0] = subx * elementsInSub[0];
				for (int z = offset[2]; z < offset[2] + elementsInSub[2]; z++) {
					for (int y = offset[1]; y < offset[1] + elementsInSub[1]; y++) {
						for (int x = offset[0]; x < offset[0] + elementsInSub[0]; x++) {

							for (int i = 0; i < 3; i++) {
								for (int j = 0; j < 3; j++) {
									for (int k = 0; k < 3; k++) {
										indices[9 * i + 3 * j + k] = nnx * nny * (2 * z + i) + nnx * (2 * y + j) + 2 * x + k;
									}
								}
							}

							tetra[0] = indices[2];
							tetra[1] = indices[6];
							tetra[2] = indices[0];
							tetra[4] = indices[20];

							tetra[8] = indices[4];
							tetra[9] = indices[3];
							tetra[11] = indices[1];
							tetra[16] = indices[11];
							tetra[17] = indices[13];
							tetra[18] = indices[10];
							mesh.pushElement(new Tetrahedron10(tetra));

							tetra[0] = indices[6];
							tetra[1] = indices[0];
							tetra[2] = indices[20];
							tetra[4] = indices[18];

							tetra[8] = indices[3];
							tetra[9] = indices[10];
							tetra[11] = indices[13];
							tetra[16] = indices[12];
							tetra[17] = indices[9];
							tetra[18] = indices[19];
							mesh.pushElement(new Tetrahedron10(tetra));

							tetra[0] = indices[24];
							tetra[1] = indices[6];
							tetra[2] = indices[20];
							tetra[4] = indices[18];

							tetra[8] = indices[15];
							tetra[9] = indices[13];
							tetra[11] = indices[22];
							tetra[16] = indices[21];
							tetra[17] = indices[12];
							tetra[18] = indices[19];
							mesh.pushElement(new Tetrahedron10(tetra));

							tetra[0] = indices[6];
							tetra[1] = indices[26];
							tetra[2] = indices[24];
							tetra[4] = indices[20];

							tetra[8] = indices[16];
							tetra[9] = indices[25];
							tetra[11] = indices[15];
							tetra[16] = indices[13];
							tetra[17] = indices[23];
							tetra[18] = indices[22];
							mesh.pushElement(new Tetrahedron10(tetra));

							tetra[0] = indices[8];
							tetra[1] = indices[26];
							tetra[2] = indices[6];
							tetra[4] = indices[20];

							tetra[8] = indices[17];
							tetra[9] = indices[16];
							tetra[11] = indices[7];
							tetra[16] = indices[14];
							tetra[17] = indices[23];
							tetra[18] = indices[13];
							mesh.pushElement(new Tetrahedron10(tetra));

							tetra[0] = indices[2];
							tetra[1] = indices[20];
							tetra[2] = indices[8];
							tetra[4] = indices[6];

							tetra[8] = indices[11];
							tetra[9] = indices[14];
							tetra[11] = indices[5];
							tetra[16] = indices[4];
							tetra[17] = indices[13];
							tetra[18] = indices[7];
							mesh.pushElement(new Tetrahedron10(tetra));
						}
					}
				}
				mesh.endPartition();
			}
		}
	}
}


void Permoncube::dirichlet(
		std::map < int, double >  & dirichlet_x,
		std::map < int, double >  & dirichlet_y,
		std::map < int, double >  & dirichlet_z,
		int *subdomains,
		int *elementsInSub)
{

	int nnx = subdomains[0] * elementsInSub[0] + 1;
	int nny = subdomains[1] * elementsInSub[1] + 1;
	int nnz = subdomains[2] * elementsInSub[2] + 1;

	idx_t index = 0;
	for (int z = 0; z < nnz; z++) {
		for (int y = 0; y < nny; y++) {
			for (int x = 0; x < nnx; x++) {
				if (z == 0){
					dirichlet_z[index] = 0.0;
				}
				if (y == 0){
					dirichlet_y[index] = 0.0;
				}
				if (x == 0){
					dirichlet_x[index] = 0.0;
				}
				index++;
			}
		}
	}
}

void dirichletTetra10(
		std::map < int, double >  & dirichlet_x,
		std::map < int, double >  & dirichlet_y,
		std::map < int, double >  & dirichlet_z,
		int *subdomains,
		int *elementsInSub)
{
	int nnx = 2 * subdomains[0] * elementsInSub[0] + 1;
	int nny = 2 * subdomains[1] * elementsInSub[1] + 1;
	int nnz = 2 * subdomains[2] * elementsInSub[2] + 1;

	for (int i = 0; i < nnx * nny; i++) {
		dirichlet_z[i] = 0.0;
		dirichlet_y[i] = 0.0;
		dirichlet_x[i] = 0.0;
	}
}

