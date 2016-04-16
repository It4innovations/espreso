#include "elementbuilder.h"

using namespace espreso::input;

ElementBuilder::ElementBuilder() {
}

ElementBuilder::~ElementBuilder() {

}

ParseError* ElementBuilder::createElement(
		std::vector<Element*> &elements) {

	std::set< eslocal > coordinates;
	int numberOfSquares = 0;
	eslocal indicies[8];
	eslocal params[Element::PARAMS_SIZE];

	for (std::list<std::pair<Face*, bool> >::iterator it =
			selectedFaces.begin(); it != selectedFaces.end(); ++it) {
		Face *face = it->first;
		coordinates.insert(face->p[0]);
		coordinates.insert(face->p[1]);
		coordinates.insert(face->p[2]);
		if (face->numberOfPoints == 4) {
			coordinates.insert(face->p[3]);
			numberOfSquares++;
		}
	}
	if (coordinates.size() == 4) {
		//Tetrahedron4
		if (selectedFaces.size() != 4) {
			std::stringstream ss;
			ss << "Element with 4 unique coordinates can not have "
					<< selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares > 0) {
			std::stringstream ss;
			ss
					<< "Element with 4 unique coordinates can not have a face with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}

		std::pair<Face*, bool> &firstFace = selectedFaces.front();
		selectedFaces.pop_front();

		coordinates.erase(firstFace.first->p[0]);
		coordinates.erase(firstFace.first->p[1]);
		coordinates.erase(firstFace.first->p[2]);

		indicies[0] = firstFace.first->p[1];
		indicies[1] = firstFace.first->p[0];
		indicies[2] = firstFace.first->p[2];
		indicies[3] = indicies[2];
		indicies[4] = *(coordinates.begin());

		Element* element = new Tetrahedron4(indicies, 8, params);
		if (firstFace.second) {
			firstFace.first->setFaceIndex(element, 0);
		}

		for (std::list<std::pair<Face*, bool> >::iterator it =
				selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it).second) {
				Face *face = (*it).first;
				if (face->containsLine(indicies[0], indicies[1]))
					face->setFaceIndex(element, 1);
				else if (face->containsLine(indicies[1], indicies[2]))
					face->setFaceIndex(element, 2);
				else
					face->setFaceIndex(element, 3);
			}
		}
		elements.push_back(element);
	} else if (coordinates.size() == 5) {
		if (selectedFaces.size() != 5) {
			std::stringstream ss;
			ss << "Element with 5 unique coordinates can not have "
					<< selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares != 1) {
			std::stringstream ss;
			ss
					<< "Element with 5 unique coordinates must have 1 face with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		std::pair<Face*, bool> firstFace;

		for (std::list<std::pair<Face*, bool> >::iterator it =
				selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it).first->numberOfPoints == 4) {
				firstFace = *it;
				break;
			}
		}

		coordinates.erase(firstFace.first->p[0]);
		coordinates.erase(firstFace.first->p[1]);
		coordinates.erase(firstFace.first->p[2]);
		coordinates.erase(firstFace.first->p[3]);

		indicies[0] = firstFace.first->p[0];
		indicies[1] = firstFace.first->p[3];
		indicies[2] = firstFace.first->p[2];
		indicies[3] = firstFace.first->p[1];
		indicies[4] = *(coordinates.begin());
		Element *element = new Pyramid5(indicies, 8, params);
		if (firstFace.second) {
			firstFace.first->setFaceIndex(element, 0);
		}

		for (std::list<std::pair<Face*, bool> >::iterator it =
				selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it).second) {
				Face *face = (*it).first;
				if (face->containsLine(indicies[0], indicies[1]))
					face->setFaceIndex(element, 1);
				else if (face->containsLine(indicies[1], indicies[2]))
					face->setFaceIndex(element, 2);
				else if (face->containsLine(indicies[2], indicies[3]))
					face->setFaceIndex(element, 3);
				else
					face->setFaceIndex(element, 4);
			}
		}
		elements.push_back(element);

	} else if (coordinates.size() == 6) {
		if (selectedFaces.size() != 5) {
			std::stringstream ss;
			ss << "Element with 6 unique coordinates can not have "
					<< selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares != 3) {
			std::stringstream ss;
			ss
					<< "Element with 6 unique coordinates must have 3 faces with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}

		std::pair<Face*, bool> firstFace;
		std::list<std::pair<Face*, bool> >::iterator firstit;
		std::pair<Face*, bool> lastFace;
		std::list<std::pair<Face*, bool> >::iterator lastit;
		bool first = true;
		for (std::list<std::pair<Face*, bool> >::iterator it =
				selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it).first->numberOfPoints == 3) {
				if (first) {
					firstFace = *it;
					firstit = it;
					first = false;
				} else {
					lastFace = *it;
					lastit = it;
				}
			}
		}
		selectedFaces.erase(firstit);
		selectedFaces.erase(lastit);

		indicies[0] = firstFace.first->p[1];
		indicies[1] = firstFace.first->p[0];
		indicies[2] = firstFace.first->p[2];
		indicies[3] = indicies[2];
		PARSE_GUARD(nextPoint(indicies[3], indicies[0], indicies[4]));
		PARSE_GUARD(nextPoint(indicies[0], indicies[1], indicies[5]));
		PARSE_GUARD(nextPoint(indicies[1], indicies[2], indicies[6]));
		indicies[7] = indicies[6];
		Element *element = new Prisma6(indicies, 8, params);
		if (firstFace.second) {
			firstFace.first->setFaceIndex(element, 3);
		}
		if (lastFace.second) {
			lastFace.first->setFaceIndex(element, 4);
		}

		for (std::list<std::pair<Face*, bool> >::iterator it =
				selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it).second) {
				Face *face = (*it).first;
				if (face->containsLine(indicies[0], indicies[1]))
					face->setFaceIndex(element, 0);
				else if (face->containsLine(indicies[1], indicies[2]))
					face->setFaceIndex(element, 1);
				else
					face->setFaceIndex(element, 2);
			}
		}
		elements.push_back(element);

	} else if (coordinates.size() == 8) {
		if (selectedFaces.size() != 6) {
			std::stringstream ss;
			ss << "Element with 8 unique coordinates can not have "
					<< selectedFaces.size() << " faces.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		if (numberOfSquares != 6) {
			std::stringstream ss;
			ss
					<< "Element with 4 unique coordinates supports only faces with 4 points.";
			return new ParseError(ss.str(), "ElementBuilder");
		}
		std::pair<Face*, bool> firstFace = selectedFaces.front();
		selectedFaces.pop_front();

		indicies[0] = firstFace.first->p[0];
		indicies[1] = firstFace.first->p[3];
		indicies[2] = firstFace.first->p[2];
		indicies[3] = firstFace.first->p[1];

		PARSE_GUARD(nextPoint(indicies[3], indicies[0], indicies[4]));
		PARSE_GUARD(nextPoint(indicies[0], indicies[1], indicies[5]));
		PARSE_GUARD(nextPoint(indicies[1], indicies[2], indicies[6]));
		PARSE_GUARD(nextPoint(indicies[2], indicies[3], indicies[7]));
		Element *element = new Hexahedron8(indicies, 8, params);

		if (firstFace.second) {
			firstFace.first->setFaceIndex(element, 4);
		}

		for (std::list<std::pair<Face*, bool> >::iterator it =
				selectedFaces.begin(); it != selectedFaces.end(); ++it) {
			if ((*it).second) {
				Face *face = (*it).first;
				if (face->containsLine(indicies[0], indicies[1]))
					face->setFaceIndex(element, 0);
				else if (face->containsLine(indicies[1], indicies[2]))
					face->setFaceIndex(element, 1);
				else if (face->containsLine(indicies[2], indicies[3]))
					face->setFaceIndex(element, 2);
				else if (face->containsLine(indicies[3], indicies[0]))
					face->setFaceIndex(element, 3);
				else
					face->setFaceIndex(element, 5);
			}
		}
		elements.push_back(element);

	} else {
		std::stringstream ss;
		ss << "Element with " << coordinates.size()
				<< " coordinates is not supported.";
		return new ParseError(ss.str(), "ElementBuilder");
	}
	return NULL;
}

ParseError* ElementBuilder::nextPoint(eslocal x, eslocal y,
eslocal &nextPoint) {
	for (std::list<std::pair<Face*, bool> >::iterator it =
			selectedFaces.begin(); it != selectedFaces.end(); ++it) {
		if ((*it).first->containsLine(x, y)) {
			return (*it).first->nextPoint(x, y, nextPoint);
		}
	}

	std::stringstream ss;
	ss << "No next point for line (" << x << "," << y << ") in Element: "
			<< this;
	return new ParseError(ss.str(), "ElementBuilder");
}
