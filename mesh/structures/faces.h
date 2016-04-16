/*
 * faces.h
 *
 *  Created on: Feb 4, 2016
 *      Author: beh01
 */

#ifndef MESH_STRUCTURES_FACES_H_
#define MESH_STRUCTURES_FACES_H_

#include <utility>
#include <vector>
#include "../elements/element.h"

namespace espreso {

/** @brief Face index contains an element and its face index */
typedef std::pair<Element*, unsigned char> FaceIndex;

class Faces {
public:
	Faces();
	virtual ~Faces();

	void reserve(std::size_t size) {
		_faces.reserve(size);
	}
	void push_back(FaceIndex &index) {
		_faces.push_back(index);
	}
	const std::vector<FaceIndex>& faceIndexes() const {
		return _faces;
	}

	const std::vector<eslocal> operator[](esglobal index) const {
		const FaceIndex &ind = _faces[index];
		return (ind.first)->getFace(ind.second);
	}
	std::size_t size() const{
		return _faces.size();
	}

private:
	/** @brief All faces indexes. */
	std::vector<FaceIndex> _faces;
};

}

#endif /* MESH_STRUCTURES_FACES_H_ */
