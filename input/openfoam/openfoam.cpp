#include "openfoam.h"
#include "foam/dictionary.h"
#include "foam/face.h"
#include "foam/elementbuilder.h"

using namespace esinput;

<<<<<<< 8627e8dd523255e077a3666ff27b8ede8aeefc14
OpenFOAM::OpenFOAM(const Options &options, int rank, int size)
{
	_projectPath = options.path;
	solveParseError(computePolyMeshPath(rank, size));
	_rank = rank;
	_size = size;
	std::cout << "Rank " << rank << ": PolyMesh path: " << _polyMeshPath
			<< std::endl;
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
	//TODO:
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
	Faces faces;
	solveParseError(parse(facesFile.getTokenizer(), faces));
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

	for (int i = 0; i < maximum; i++) {
		elementBuilders.push_back(new ElementBuilder(&faces));
	}

	esglobal face = 0;
	for (std::vector<esglobal>::iterator it = owner.begin(); it != owner.end();
			++it) {
		elementBuilders[*it]->add(face);
		face++;
	}

	FoamFile neighbourFile(_polyMeshPath + "neighbour");
	std::vector< esglobal > neighbour;

	solveParseError(parse(neighbourFile.getTokenizer(), neighbour));

	face = 0;
	for (std::vector<esglobal>::iterator it = neighbour.begin();
			it != neighbour.end(); ++it) {
		elementBuilders[*it]->add(face);
		face++;
	}
	for (std::vector<ElementBuilder*>::iterator it = elementBuilders.begin();
			it != elementBuilders.end(); ++it) {
		solveParseError((*it)->createElement(elements));
		delete *it;
	}

}

void OpenFOAM::boundaryConditions(mesh::Coordinates &coordinates) {

}

void OpenFOAM::clusterBoundaries(mesh::Mesh &mesh,
		mesh::Boundaries &boundaries) {

}

