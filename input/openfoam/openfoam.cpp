#include "openfoam.h"

using namespace esinput;

OpenFOAM::OpenFOAM(const Options &options, int rank, int size)
{
	_projectPath = options.path;
	solveParseError(computePolyMeshPath(rank, size));
	_rank = rank;
	_size = size;
	std::stringstream ss;
	ss<< "Rank " << rank << "- PolyMesh path: " << _polyMeshPath
			<< std::endl;
	std::cout << ss.str();
}

void OpenFOAM::solveParseError(ParseError *error) {
	if (error != NULL) {
		error->print();
		delete error;
		exit(EXIT_FAILURE);
	}
}

ParseError* OpenFOAM::computePolyMeshPath(int rank, int size) {
	if (size > 1) {

		std::string decomposePar = _projectPath + "/system/decomposeParDict";
		std::string mesh;

		if (FoamFile::exists(decomposePar)) {
			FoamFile file(decomposePar);
			Dictionary dictionary;
			PARSE_GUARD(dictionary.parseTopLevel(file.getTokenizer()));
			int numberOfSubdomains;
			PARSE_GUARD(
					dictionary.readEntry("numberOfSubdomains",
							numberOfSubdomains));
			if (numberOfSubdomains != size) {
				std::stringstream ss;
				ss << "Task is decopmosed to " << numberOfSubdomains
						<< " sub-domains. But there is " << size
						<< " processes.";
				return new ParseError(ss.str(), "Loader");
			}
			std::stringstream ss;
			ss << _projectPath << "/processor" << rank << "/constant/polyMesh/";
			_polyMeshPath = ss.str();
			return NULL;
		}
		std::stringstream ss;
		ss << "There is no file: " << decomposePar << "; but there are " << size
				<< " processes.";
		return new ParseError(ss.str(), "Loader");
	}
	_polyMeshPath = _projectPath + "/constant/polyMesh/";
	return NULL;
}

void OpenFOAM::points(mesh::Coordinates &coordinates) {
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
			std::stringstream ss;
			ss << "The number of points(" << points.size()
					<< ") is not the same as the number of points(";
			ss << addressing.size() << ") in the file:" << _polyMeshPath
					<< "pointProcAddressing\n";
			std::cerr << ss.str();
			exit(EXIT_FAILURE);
		}
		std::vector< esglobal >::iterator it_addressing = addressing.begin();
		it_addressing = addressing.begin();
		for (std::vector<mesh::Point>::iterator it = points.begin();
				it != points.end(); ++it) {
			coordinates.add(*it, counter, *it_addressing);
			++it_addressing;
		}
	} else {
		for (std::vector<mesh::Point>::iterator it = points.begin();
				it != points.end(); ++it) {
			coordinates.add(*it, counter, counter);
		}
	}
}

void OpenFOAM::elements(std::vector<mesh::Element*> &elements) {
	FoamFile facesFile(_polyMeshPath + "faces");
	solveParseError(parse(facesFile.getTokenizer(), _faces));
	FoamFile ownerFile(_polyMeshPath + "owner");
	std::vector< esglobal > owner;
	solveParseError(parse(ownerFile.getTokenizer(), owner));

	esglobal maximum = 0;
	for (std::vector<esglobal>::iterator it = owner.begin(); it != owner.end();
			++it) {
		if (*it + 1 > maximum) {
			maximum = *it + 1;
		}
	}

	std::vector<ElementBuilder*> elementBuilders;
	elements.reserve(maximum);
	elementBuilders.reserve(maximum);

	//#pragma omp for
	for (int i = 0; i < maximum; i++) {
		elementBuilders.push_back(new ElementBuilder());
	}

	esglobal face = 0;
	for (std::vector<esglobal>::iterator it = owner.begin(); it != owner.end();
			++it) {
		elementBuilders[*it]->add(&(_faces.at(face)), true);
		face++;
	}

	FoamFile neighbourFile(_polyMeshPath + "neighbour");
	std::vector< esglobal > neighbour;

	solveParseError(parse(neighbourFile.getTokenizer(), neighbour));

	face = 0;
	for (std::vector<esglobal>::iterator it = neighbour.begin();
			it != neighbour.end(); ++it) {
		elementBuilders[*it]->add(&(_faces.at(face)), false);
		face++;
	}
	for (std::vector<ElementBuilder*>::iterator it = elementBuilders.begin();
			it != elementBuilders.end(); ++it) {
		solveParseError((*it)->createElement(elements));
		delete *it;
	}
	//reads also cell zones
	FoamFile cellZonesFile(_polyMeshPath + "cellZones");
	solveParseError(parse(cellZonesFile.getTokenizer(), _cellZones));
	/*for (std::vector<CellZone>::iterator it = _cellZones.begin();
				it != _cellZones.end(); ++it) {
		std::cout<<*it;
	}*/
}

void OpenFOAM::faces(mesh::Faces &faces) {
	for (std::vector<Face>::iterator it = _faces.begin(); it != _faces.end();
			++it) {
		faces.push_back((*it).getFaceIndex());
	}
	_faces.clear();
}

void OpenFOAM::boundaryConditions(mesh::Coordinates &coordinates) {

}

void OpenFOAM::clusterBoundaries(mesh::Mesh &mesh,
		mesh::Boundaries &boundaries) {
	boundaries.resize(mesh.coordinates().clusterSize());
	for (size_t i = 0; i < mesh.coordinates().clusterSize(); i++) {
		boundaries[i].insert(_rank);
	}
	if (_size > 1) {

/*		for (int i = 0; i < mesh.faces().size(); i++) {
			std::cout << "Face " << i;
			std::vector<eslocal> face = mesh.faces()[i];
			std::cout << " (";
			for (std::vector<eslocal>::iterator it = face.begin();
					it != face.end(); ++it) {
				std::cout << *it << ",";
			}
			std::cout << ")\n";

		}*/

		FoamFile boundaryFile(_polyMeshPath + "boundary");
		std::vector<Dictionary> boundary;
		solveParseError(parse(boundaryFile.getTokenizer(), boundary));

		for (std::vector<Dictionary>::iterator it = boundary.begin();
				it != boundary.end(); ++it) {
			if ((*it).getName().find("procBoundary") == 0) {
				int myProcNo;
				solveParseError((*it).readEntry("myProcNo", myProcNo));
				if (myProcNo != _rank) {
					std::stringstream ss;
					ss << "Boundary for rank: " << myProcNo
							<< ", but opened in process: " << _rank;
					std::cerr << ss.str();
					exit(EXIT_FAILURE);
				}

				int neighbProcNo;
				solveParseError((*it).readEntry("neighbProcNo", neighbProcNo));
				int nFaces;
				solveParseError((*it).readEntry("nFaces", nFaces));
				int startFace;
				solveParseError((*it).readEntry("startFace", startFace));
				for (int i = 0; i < nFaces; i++) {
					std::vector<eslocal> face = mesh.faces()[startFace + i];
				//	std::cout << "Face: ";
					for (std::vector<eslocal>::iterator it = face.begin();
							it != face.end(); ++it) {
				//		std::cout << *it << ",";
						boundaries[*it].insert(neighbProcNo);
					}
				//	std::cout << "\n";
				}
				/*std::cout << "Process: " << myProcNo << " " << neighbProcNo
						<< " " << startFace << " " << nFaces << "\n";
				*/
			}
		}
	}
}

