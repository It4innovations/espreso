
#include "openfoam.h"
#include "foam/dictionary.h"

using namespace esinput;

OpenFOAM::OpenFOAM(const Options &options, int rank, int size)
{
	_projectPath = options.path;
	solveParseError(computePolyMeshPath(rank, size));
	_rank = rank;
	_size = size;
	std::cout<<"Rank "<<rank<<": PolyMesh path: "<<_polyMeshPath<<std::endl;
}

void OpenFOAM::solveParseError(ParseError *error) {
	if (error!=NULL) {
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
			PARSE_GUARD(dictionary.readEntry("numberOfSubdomains", numberOfSubdomains));
			if (numberOfSubdomains!=size) {
				std::stringstream ss;
				ss<<"Task is decopmosed to "<<numberOfSubdomains<<" sub-domains. But there is "<<size<<" processes.";
				return new ParseError(ss.str(), "Loader");
			}
			std::stringstream ss;
			ss<<_projectPath<<"/processor"<<rank<<"/constant/polyMesh/";
			_polyMeshPath=ss.str();
			return NULL;
    	}
		std::stringstream ss;
		ss<<"There is no file: "<<decomposePar<<"; but there are "<<size<<" processes.";
		return new ParseError(ss.str(), "Loader");
	}
    _polyMeshPath = _projectPath + "/constant/polyMesh/";
	return NULL;
}

void OpenFOAM::points(mesh::Coordinates &coordinates)
{
	FoamFile pointsFile(_polyMeshPath+"points");
	//TODO:
	Points points;
	solveParseError(parse(pointsFile.getTokenizer(), points));
	coordinates.reserve(points.size());
	eslocal counter =0;
	if (_size>1) {
		FoamFile addressingFile(_polyMeshPath+"pointProcAddressing");
		std::vector< esglobal > addressing;
		solveParseError(parse(addressingFile.getTokenizer(), addressing));
		std::vector< esglobal >::iterator it_addressing = addressing.begin();
		it_addressing = addressing.begin();
		for (std::vector< mesh::Point >::iterator it = points.begin() ; it != points.end(); ++it) {
			coordinates.add(*it, counter, *it_addressing);
			++it_addressing;
		}
	}else {
		for (std::vector< mesh::Point >::iterator it = points.begin() ; it != points.end(); ++it) {
			coordinates.add(*it, counter, counter);
		}
	}
}


void OpenFOAM::elements(std::vector<mesh::Element*> &elements)
{


}

void OpenFOAM::boundaryConditions(mesh::Coordinates &coordinates)
{

}

void OpenFOAM::clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries)
{

}




