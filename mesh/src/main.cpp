#include "esmesh.h"

int main(int argc, char** argv)
{
	int partsCount = 4;
	int fixPointsCount = 8;

	Coordinates coords("matrices/HEX/5/coord");
	Mesh mesh(coords);
	mesh = Mesh("matrices/HEX/5/elem", coords, partsCount, fixPointsCount);
	Boundaries boundaries(mesh, coords);
	Faces faces(mesh, coords);
	Corners corners(faces.getFaces(), coords);

	mesh.saveVTK();
}
