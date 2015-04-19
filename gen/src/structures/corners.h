#ifndef CORNERS_H_
#define CORNERS_H_

#include "../elements/elements.h"
#include "mesh.h"

class Corners
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Corners &c);

	Corners(const BoundaryFaces &faces, const Coordinates &coordinates);
	~Corners();

private:
	/** @brief Reference to coordinates. */
	const Coordinates &_coordinates;

	/** @brief Array that stores all corners elements. */
	BoundaryLines _lines;
};



#endif /* CORNERS_H_ */
