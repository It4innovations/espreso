
#include "openfoam.h"
#include "foam/dictionary.h"

using namespace esinput;

OpenFOAM::OpenFOAM(const Options &options, int rank, int size)
{
	_projectPath = options.path;
	ParseError *error = computePolyMeshPath(rank, size);
	if (error!=NULL) {
		error->print();
		delete error;
		exit(EXIT_FAILURE);
	}
	std::cout<<"Rank "<<rank<<": PolyMesh path: "<<_polyMeshPath;
}

ParseError* OpenFOAM::computePolyMeshPath(int rank, int size) {
	std::string decomposePar = _projectPath + "/system/decomposeParDict";
	std::string mesh;
	if (FoamFile::exists(decomposePar)) {
		FoamFile file(decomposePar);
	    Dictionary dictionary;
	    PARSE_GUARD(dictionary.parseTopLevel(file.getTokenizer()));
	    int numberOfSubdomains;
        PARSE_GUARD(dictionary.readEntry("numberOfSubdomains", numberOfSubdomains));
        std::stringstream ss;
        ss<<_projectPath<<"/processor"<<rank<<"/constant/polyMesh/";
        _polyMeshPath=ss.str();
    }else{
        _polyMeshPath = _projectPath + "/constant/polyMesh/";
    }
	return NULL;
}

void OpenFOAM::points(mesh::Coordinates &coordinates)
{
	FoamFile pointsFile(mesh+"points");
	Points points;
	parse(pointsFile.getTokenizer(), points);
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




