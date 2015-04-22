#include "esmesh.h"

int main(int argc, char** argv)
{
	int partsCount = 4;
	int fixPointsCount = 8;

	Coordinates coords("matrices/HEX/15/coord");
	Mesh mesh("matrices/HEX/15/elem", coords, partsCount, fixPointsCount);
	Boundaries boundaries(mesh, coords);
	Faces faces(mesh, coords);
	Corners corners(faces.getFaces(), coords);
}
