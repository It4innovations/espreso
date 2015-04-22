#include "corners.h"

Corners::Corners(const BoundaryFaces &faces, const Coordinates &coordinates): _coordinates(coordinates)
{
	BoundaryFaces::const_iterator it;
	int parts[2];
	for (it = faces.begin(); it != faces.end(); ++it) {
		if (it->second.first != it->second.second) {
			parts[0] = it->second.first;
			parts[1] = it->second.second;
			it->first->fillLines(_lines, parts);
		}
	}
}

Corners::~Corners()
{
	BoundaryLines::iterator it;
	for (it = _lines.begin(); it != _lines.end(); ++it) {
		delete it->first;
	}
}

std::ostream& operator<<(std::ostream& os, const Corners &c)
{
	size_t inner = 0, border = 0;
	BoundaryLines::const_iterator it;
	std::set<int>::const_iterator sit;
	for (it = c._lines.begin(); it != c._lines.end(); ++it) {
		if (it->second.size() > 2) {
			os << *(it->first) << ": ";
			for (sit = it->second.begin(); sit != it->second.end(); ++sit) {
				os << *sit << " ";
			}
			os << "\n";
			border++;
		} else {
			inner++;
		}
	}
	os << "inner: " << inner << "\n";
	os << "border: " << border << "\n";
	return os;
}



