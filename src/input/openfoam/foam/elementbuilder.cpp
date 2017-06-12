
#include "elementbuilder.h"
#include "../../../mesh/structures/coordinates.h"
#include "../../../mesh/elements/volume/hexahedron8.h"
#include "../../../mesh/elements/volume/hexahedron8.h"
#include "../../../mesh/elements/volume/tetrahedron4.h"
#include "../../../mesh/elements/volume/prisma6.h"
#include "../../../mesh/elements/volume/pyramid5.h"

using namespace espreso::input;

ElementBuilder::ElementBuilder() {
}

ElementBuilder::~ElementBuilder() {

}

ParseError* ElementBuilder::createElement(Element *&element, const Coordinates &coordinates) {

	std::set< eslocal > nodes;
	int numberOfSquares = 0;
	eslocal indicies[8];
	std::vector<eslocal> params(Element::PARAMS_SIZE);

	for (std::list<Face* >::iterator  it = selectedFaces.begin(); it != selectedFaces.end(); ++it) {
		Face *face = *it;
		nodes.insert(face->p[0]);
		nodes.insert(face->p[1]);
		nodes.insert(face->p[2]);
		if (face->numberOfPoints == 4) {
			nodes.insert(face->p[3]);
			numberOfSquares++;
		}
	}
	if (nodes.size() == 4) {
		//Tetrahedron4
		if (selectedFaces.size() != 4) {
			std::stringstream ss;
			ss << "Element with 4 unique coordinates can not have " << selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares > 0) {
			std::stringstream ss;
			ss << "Element with 4 unique coordinates can not have a face with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}

		Face* firstFace = selectedFaces.front();
		selectedFaces.pop_front();

		nodes.erase(firstFace->p[0]);
		nodes.erase(firstFace->p[1]);
		nodes.erase(firstFace->p[2]);

		indicies[0] = firstFace->p[1];
		indicies[1] = firstFace->p[2];
		indicies[2] = firstFace->p[0];
		indicies[3] = *(nodes.begin());

		Point normal = Point::cross(coordinates[indicies[1]] - coordinates[indicies[0]], coordinates[indicies[2]] - coordinates[indicies[0]]);
		if (normal * (coordinates[indicies[3]] - coordinates[indicies[0]]) < 0) {
			std::swap(indicies[1], indicies[2]);
		}

		element = new Tetrahedron4(indicies, 4, params.data());

	} else if (nodes.size() == 5) {
		if (selectedFaces.size() != 5) {
			std::stringstream ss;
			ss << "Element with 5 unique nodes can not have " << selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares != 1) {
			std::stringstream ss;
			ss << "Element with 5 unique coordinates must have 1 face with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		Face* firstFace = NULL;

		for (std::list<Face* >::iterator it = selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it)->numberOfPoints == 4) {
				firstFace = *it;
				break;
			}
		}

		nodes.erase(firstFace->p[0]);
		nodes.erase(firstFace->p[1]);
		nodes.erase(firstFace->p[2]);
		nodes.erase(firstFace->p[3]);

		indicies[0] = firstFace->p[0];
		indicies[1] = firstFace->p[3];
		indicies[2] = firstFace->p[2];
		indicies[3] = firstFace->p[1];
		indicies[4] = *(nodes.begin());

		Point normal = Point::cross(coordinates[indicies[1]] - coordinates[indicies[0]], coordinates[indicies[2]] - coordinates[indicies[0]]);
		if (normal * (coordinates[indicies[4]] - coordinates[indicies[0]]) < 0) {
			std::swap(indicies[1], indicies[3]);
		}

		element = new Pyramid5(indicies, 5, params.data());

	} else if (nodes.size() == 6) {
		if (selectedFaces.size() != 5) {
			std::stringstream ss;
			ss << "Element with 6 unique nodes can not have " << selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares != 3) {
			std::stringstream ss;
			ss << "Element with 6 unique coordinates must have 3 faces with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}

		Face* firstFace = NULL;
		std::list< Face* >::iterator firstit;
		std::list< Face* >::iterator lastit;
		bool first = true;
		for (std::list<Face* >::iterator it = selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it)->numberOfPoints == 3) {
				if (first) {
					firstFace = *it;
					firstit = it;
					first = false;
				} else {
					lastit = it;
				}
			}
		}
		selectedFaces.erase(firstit);
		selectedFaces.erase(lastit);

		indicies[0] = firstFace->p[1];
		indicies[1] = firstFace->p[0];
		indicies[2] = firstFace->p[2];
		PARSE_GUARD(nextPoint(indicies[2], indicies[0], indicies[3]));
		PARSE_GUARD(nextPoint(indicies[0], indicies[1], indicies[4]));
		PARSE_GUARD(nextPoint(indicies[1], indicies[2], indicies[5]));

		Point normal = Point::cross(coordinates[indicies[1]] - coordinates[indicies[0]], coordinates[indicies[2]] - coordinates[indicies[0]]);
		if (normal * (coordinates[indicies[3]] - coordinates[indicies[0]]) < 0) {
			std::swap(indicies[0], indicies[3]);
			std::swap(indicies[1], indicies[4]);
			std::swap(indicies[2], indicies[5]);
		}

		element = new Prisma6(indicies, 6, params.data());

	} else if (nodes.size() == 8) {
		if (selectedFaces.size() != 6) {
			std::stringstream ss;
			ss << "Element with 8 unique nodes can not have " << selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares != 6) {
			std::stringstream ss;
			ss << "Element with 4 unique coordinates supports only faces with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		Face* firstFace = selectedFaces.front();
		selectedFaces.pop_front();

		indicies[0] = firstFace->p[0];
		indicies[1] = firstFace->p[3];
		indicies[2] = firstFace->p[2];
		indicies[3] = firstFace->p[1];

		PARSE_GUARD(nextPoint(indicies[3], indicies[0], indicies[4]));
		PARSE_GUARD(nextPoint(indicies[0], indicies[1], indicies[5]));
		PARSE_GUARD(nextPoint(indicies[1], indicies[2], indicies[6]));
		PARSE_GUARD(nextPoint(indicies[2], indicies[3], indicies[7]));

		Point center, centerup, centerdown;
		for (size_t i = 0; i < 4; i++) {
			centerup += coordinates[indicies[i]];
			centerdown += coordinates[indicies[i + 4]];
		}
		center = (centerup + centerdown) / 8.0;
		centerdown /= 4.0;
		centerup /= 4.0;

		Point normala = Point::cross(coordinates[indicies[1]] - coordinates[indicies[0]], coordinates[indicies[2]] - coordinates[indicies[0]]);
		Point normalb = Point::cross(coordinates[indicies[2]] - coordinates[indicies[0]], coordinates[indicies[3]] - coordinates[indicies[0]]);
		if (((normala + normalb) / 2.0) * (centerup - centerdown) >= 0) {
			std::swap(indicies[0], indicies[4]);
			std::swap(indicies[1], indicies[5]);
			std::swap(indicies[2], indicies[6]);
			std::swap(indicies[3], indicies[7]);
		}

		element = new Hexahedron8(indicies, 8, params.data());
	} else {
		std::stringstream ss;
		ss << "Element with " << nodes.size() << " coordinates is not supported.";
		return new ParseError(ss.str(), "ElementBuilder");
	}
	return NULL;
}

ParseError* ElementBuilder::nextPoint(eslocal x, eslocal y, eslocal &nextPoint) {
	for (std::list<Face* >::iterator it = selectedFaces.begin(); it != selectedFaces.end(); ++it) {
		if ((*it)->containsLine(x, y)) {
			return (*it)->nextPoint(x, y, nextPoint);
		}
	}

	std::stringstream ss;
	ss << "No next point for line (" << x << "," << y << ") in Element: " << this;
	return new ParseError(ss.str(), "ElementBuilder");
}
