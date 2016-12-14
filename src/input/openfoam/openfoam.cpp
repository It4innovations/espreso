#include "openfoam.h"

using namespace espreso::input;

OpenFOAM::OpenFOAM(const GlobalConfiguration &configuration, Mesh &mesh, int rank, int size) :
		Loader(configuration, mesh), _openfoam(configuration.openfoam)
{
	_projectPath = _openfoam.path;
	solveParseError(computePolyMeshPath(rank, size));
	_rank = rank;
	_size = size;
}

void OpenFOAM::solveParseError(ParseError *error)
{
	if (error != NULL) {
		error->print();
		delete error;
		exit(EXIT_FAILURE);
	}
}

ParseError* OpenFOAM::computePolyMeshPath(int rank, int size)
{
	if (size > 1) {
		std::string decomposePar = _projectPath + "/system/decomposeParDict";
		std::string mesh;

		if (FoamFile::exists(decomposePar)) {
			FoamFile file(decomposePar);
			Dictionary dictionary;
			PARSE_GUARD(dictionary.parseTopLevel(file.getTokenizer()));
			int numberOfSubdomains;
			PARSE_GUARD(dictionary.readEntry("numberOfSubdomains", numberOfSubdomains));
			if (numberOfSubdomains != size) {
				std::stringstream ss;
				ss << "Task is decopmosed to " << numberOfSubdomains << " sub-domains. But there is " << size << " processes.";
				return new ParseError(ss.str(), "Loader");
			}
			std::stringstream ss;
			ss << _projectPath << "/processor" << rank << "/constant/polyMesh/";
			_polyMeshPath = ss.str();
			return NULL;
		}
		std::stringstream ss;
		ss << "There is no file: " << decomposePar << "; but there are " << size << " processes.";
		return new ParseError(ss.str(), "Loader");
	}
	_polyMeshPath = _projectPath + "/constant/polyMesh/";
	return NULL;
}

void OpenFOAM::points(Coordinates &coordinates)
{
	FoamFile pointsFile(_polyMeshPath + "points");
	Points points;
	solveParseError(parse(pointsFile.getTokenizer(), points));
	coordinates.reserve(points.size());
	eslocal counter = 0;
	if (_size > 1) {
		FoamFile addressingFile(_polyMeshPath + "pointProcAddressing");
		std::vector< esglobal > addressing;
		solveParseError(parse(addressingFile.getTokenizer(), addressing));
		if (addressing.size() != points.size()) {
			ESINFO(GLOBAL_ERROR)
					<< "The number of points(" << points.size() << ") is not the same as the number of points("
					<< addressing.size() << ") in the file:" << _polyMeshPath << "pointProcAddressing\n";
		}
		std::vector< esglobal >::iterator it_addressing = addressing.begin();
		for (std::vector<Point>::iterator it = points.begin(); it != points.end(); ++it, counter++) {
			coordinates.add(*it, counter, *it_addressing);
			++it_addressing;
		}
	} else {
		for (std::vector<Point>::iterator it = points.begin(); it != points.end(); ++it, counter++) {
			coordinates.add(*it, counter, counter);
		}
	}
}

void OpenFOAM::elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	FoamFile facesFile(_polyMeshPath + "faces");
	std::vector<Face> _faces;
	solveParseError(parse(facesFile.getTokenizer(), _faces));

	for (std::vector<Face>::iterator it = _faces.begin(); it != _faces.end(); ++it) {
		Element *faceElement;
		switch ((*it).numberOfPoints) {
		case 3:
			faceElement = new Triangle3((*it).p);
			break;
		case 4:
			faceElement = new Square4((*it).p);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Faces with " << (*it).numberOfPoints << " points are not supported.";
		}
		faces.push_back(faceElement);
		(*it).faceElement = faceElement;
	}

	FoamFile ownerFile(_polyMeshPath + "owner");
	std::vector< esglobal > owner;
	solveParseError(parse(ownerFile.getTokenizer(), owner));

	esglobal maximum = 0;
	for (std::vector<esglobal>::iterator it = owner.begin(); it != owner.end(); ++it) {
		if (*it + 1 > maximum) {
			maximum = *it + 1;
		}
	}

	std::vector<ElementBuilder*> elementBuilders;
	elements.reserve(maximum);
	elementBuilders.reserve(maximum);

	for (int i = 0; i < maximum; i++) {
		elementBuilders.push_back(new ElementBuilder());
	}

	esglobal face = 0;
	for (std::vector<esglobal>::iterator it = owner.begin(); it != owner.end(); ++it) {
		elementBuilders[*it]->add(&_faces[face], true);
		face++;
	}

	FoamFile neighbourFile(_polyMeshPath + "neighbour");
	std::vector< esglobal > neighbour;

	solveParseError(parse(neighbourFile.getTokenizer(), neighbour));

	face = 0;
	for (std::vector<esglobal>::iterator it = neighbour.begin(); it != neighbour.end(); ++it) {
		elementBuilders[*it]->add(&_faces[face], false);
		face++;
	}
	for (std::vector<ElementBuilder*>::iterator it = elementBuilders.begin(); it != elementBuilders.end(); ++it) {
		VolumeElement *element = NULL;
		solveParseError((*it)->createElement(element));
		if (element == NULL) {
			ESINFO(GLOBAL_ERROR) << "Unrecognized element form face: " << (*it) << "\n";
		} else {
			elements.push_back(element);
		}
		delete *it;
	}

	for (size_t i = 0; i < owner.size(); i++) {
		elements[owner[i]]->addFace(faces[i]);
	}
	for (size_t i = 0; i < neighbour.size(); i++) {
		elements[neighbour[i]]->addFace(faces[i]);
	}
}

void OpenFOAM::materials(std::vector<Material> &materials)
{
	materials.resize(1, Material(mesh.coordinates()));
}

void OpenFOAM::regions(
		std::vector<Evaluator*> &evaluators,
		std::vector<Region> &regions,
		std::vector<Element*> &elements,
		std::vector<Element*> &faces,
		std::vector<Element*> &edges,
		std::vector<Element*> &nodes)
{
	//reads boundary
	FoamFile boundaryFile(_polyMeshPath + "boundary");
	std::vector<Dictionary> boundary;
	solveParseError(parse(boundaryFile.getTokenizer(), boundary));

	for (std::vector<Dictionary>::iterator it = boundary.begin(); it != boundary.end(); ++it) {
		eslocal nFaces = 0;
		solveParseError((*it).readEntry("nFaces", nFaces));
		eslocal startFace = 0;
		solveParseError((*it).readEntry("startFace", startFace));

		if ((*it).getName().find("procBoundary") != 0) {
			regions.push_back(Region());
			Region *region = &(regions[regions.size() - 1]);
			region->name = (*it).getName();
			region->elements.resize(nFaces);
			memcpy(&region->elements[0], &faces[startFace], nFaces * sizeof(Element*));
		}
	}
	//reads cell zones
	if (!FoamFile::checkFileType(_polyMeshPath + "cellZones").empty()) {
		FoamFile cellZonesFile(_polyMeshPath + "cellZones");
		std::vector<CellZone> _cellZones;
		solveParseError(parse(cellZonesFile.getTokenizer(), _cellZones));

		for (auto cellZone : _cellZones) {
			regions.push_back(Region());
			Region *region = &(regions[regions.size() - 1]);
			region->name = cellZone.getName();
			for (auto index : cellZone.elementIndexes()) {
				region->elements.push_back(elements[index]);
			}
		}
	}
	if (!FoamFile::checkFileType(_polyMeshPath + "faceZones").empty()) {
		//reads face zones
		FoamFile faceZonesFile(_polyMeshPath + "faceZones");
		std::vector<FaceZone> _faceZones;
		solveParseError(parse(faceZonesFile.getTokenizer(), _faceZones));
		for (auto faceZone : _faceZones) {
			regions.push_back(Region());
			Region *region = &(regions[regions.size() - 1]);
			region->name = faceZone.getName();
			for (auto index : faceZone.elementIndexes()) {
				region->elements.push_back(faces[index]);
			}
		}
	}
	if (!FoamFile::checkFileType(_polyMeshPath + "pointZones").empty()) {
		//reads point zones
		FoamFile pointZonesFile(_polyMeshPath + "pointZones");
		std::vector<PointZone> _pointZones;
		solveParseError(parse(pointZonesFile.getTokenizer(), _pointZones));
		for (auto pointZone : _pointZones) {
			regions.push_back(Region());
			Region *region = &(regions[regions.size() - 1]);
			region->name = pointZone.getName();
			for (auto index : pointZone.elementIndexes()) {
				region->elements.push_back(nodes[index]);
			}
		}
	}
}

void OpenFOAM::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	FoamFile boundaryFile(_polyMeshPath + "boundary");
	std::vector<Dictionary> boundary;
	solveParseError(parse(boundaryFile.getTokenizer(), boundary));

	int myProcNo = 0;
	int neighbProcNo = 0;
	int lastNeighbProcNo = -1;

	for (std::vector<Dictionary>::iterator it = boundary.begin(); it != boundary.end(); ++it) {
		eslocal nFaces = 0;
		solveParseError((*it).readEntry("nFaces", nFaces));
		eslocal startFace = 0;
		solveParseError((*it).readEntry("startFace", startFace));

		if ((*it).getName().find("procBoundary") == 0) {
			solveParseError((*it).readEntry("myProcNo", myProcNo));
			if (myProcNo != _rank) {
				ESINFO(GLOBAL_ERROR) << "Boundary for rank: " << myProcNo << ", but opened in process: " << _rank;
			}
			solveParseError((*it).readEntry("neighbProcNo", neighbProcNo));

			if (lastNeighbProcNo < _rank && _rank < neighbProcNo) {
				for (size_t n = 0; n < nodes.size(); n++) {
					nodes[n]->clusters().push_back(_rank);
				}
			}

			for (eslocal i = startFace; i < startFace + nFaces; i++) {
				for (size_t n = 0; n < faces[i]->nodes(); n++) {
					if (!nodes[faces[i]->node(n)]->clusters().size() || nodes[faces[i]->node(n)]->clusters().back() != neighbProcNo) {
						nodes[faces[i]->node(n)]->clusters().push_back(neighbProcNo);
					}
				}
			}

			lastNeighbProcNo = neighbProcNo;
			neighbours.push_back(neighbProcNo);
		}
	}

	if (lastNeighbProcNo < _rank) {
		for (size_t n = 0; n < nodes.size(); n++) {
			nodes[n]->clusters().push_back(_rank);
		}
	}
}

