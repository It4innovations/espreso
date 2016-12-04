#ifndef CELLZONE_H
#define CELLZONE_H

#include "dictionary.h"
#include <vector>
#include <string>
#include "../../loader.h"

namespace espreso {
namespace input {

class Zone {
public:
	Zone(std::string containerIndexesName);

	const std::string& getName() const {
		return name;
	}

	ParseError* loadFromDictionary(Dictionary &dictionary);

	const std::vector<eslocal >& elementIndexes() const {
		return _elementIndexes;
	}

	friend inline std::ostream& operator<<(std::ostream& os,
			const Zone& obj) {
		// write obj to stream
		os << obj.getName() << "{\n";
		os << "cellLabels " << obj._elementIndexes.size() << "(\n";
		for (std::vector< eslocal >::const_iterator it = obj._elementIndexes.begin();
				it != obj._elementIndexes.end(); ++it) {
			os << *it << "\n";
		}
		os << ")};\n";
		return os;
	}

protected:
	std::string name;
	std::string containerIndexesName;
	std::vector< eslocal > _elementIndexes;
};

class CellZone : public Zone {
public:
	CellZone() : Zone("cellLabels") {

	}

};

class FaceZone : public Zone {
public:
	FaceZone() : Zone("faceLabels") {

	}

};

ParseError* parse(Tokenizer &ts, CellZone &cellZone);
ParseError* parse(Tokenizer &ts, FaceZone &faceZone);
ParseError* parse(Tokenizer &ts, Zone &Zone);

}
}

#endif // CELLZONE_H
