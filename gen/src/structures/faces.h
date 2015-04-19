#ifndef FACES_H_
#define FACES_H_

#include "../elements/elements.h"
#include "mesh.h"

class Faces
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Faces &f);

	Faces(const Mesh &mesh, const Coordinates &coordinates);
	Faces(const Mesh &mesh, const Coordinates &coordinates, const BoundaryNodes &nodes);
	~Faces();

	const BoundaryFaces& getFaces() const
	{
		return _faces;
	}

private:
	/** @brief Reference to coordinates. */
	const Coordinates &_coordinates;

	/** @brief Array that stores all elements of the mesh. */
	BoundaryFaces _faces;

	/** @brief Keeps mapping of nodes to mesh parts. */
	BoundaryNodes _boundaries;
};



#endif /* FACES_H_ */
