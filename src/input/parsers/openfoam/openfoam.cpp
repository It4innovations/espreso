
#include "openfoam.h"

#include "parser/faceList.h"
#include "parser/labelList.h"
#include "parser/vectorField.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "config/ecf/input/input.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "input/builders/builder.h"

using namespace espreso;

InputOpenFoam::InputOpenFoam()
{
	FoamFile::init();
}

InputOpenFoam::~InputOpenFoam()
{
	FoamFile::finish();
}

void InputOpenFoam::load(const InputConfiguration &configuration)
{
	eslog::startln("OPENFOAM PARSER: STARTED", "OPENFOAM PARSER");

	InputFilePack inputPack;
	AsyncFilePack asyncPack;

	OpenFOAMVectorField points(inputPack.add(configuration.path + "/constant/polyMesh/points"));
	OpenFOAMFaceList faces(inputPack.add(configuration.path + "/constant/polyMesh/faces"));

	OpenFOAMLabelList owner(asyncPack.add(configuration.path + "/constant/polyMesh/owner"));
	OpenFOAMLabelList neighbour(asyncPack.add(configuration.path + "/constant/polyMesh/neighbour"));

	inputPack.prepare();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY READER PREPARED");

	inputPack.read();
	profiler::synccheckpoint("geometry_read");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY READ");

	points.scan(); faces.scan();
	FoamFile::synchronize({ &points, &faces });
	profiler::synccheckpoint("scan");

	// compute better distribution for binary data according to scanned parameters
	asyncPack.setTotalSizes();
	owner.input->setDistribution(tarray<size_t>::distribute(info::mpi::size, owner.input->totalSize));
	neighbour.input->setDistribution(tarray<size_t>::distribute(info::mpi::size, neighbour.input->totalSize));
//	asyncPack.iread([&] () {
//		owner.scan(); neighbour.scan();
//		FoamFile::synchronize({ &owner, &neighbour });
//		owner.parse(mesh.elements->owner); neighbour.parse(mesh.elements->neighbor);
//		std::vector<std::string> info = Parser::split(owner.header.note, " ");
//		for (size_t i = 0; i < info.size(); ++i) {
//			if (StringCompare::caseSensitivePreffix("nCells", info[i])) {
//				mesh.elements->elements = std::strtol(Parser::split(info[i], ":")[1].c_str(), nullptr, 10);
//			}
//		}
//	});

	points.parse(mesh.nodes->coordinates);
	faces.parse(mesh.elements->etype, mesh.elements->enodes);

	asyncPack.wait();
	profiler::synccheckpoint("parse");
	eslog::checkpointln("OPENFOAM PARSER: GEOMETRY PARSED");

	eslog::endln("OPENFOAM PARSER: PARSED");
}

void InputOpenFoam::build(Mesh &mesh)
{
	builder::buildOrderedFVM(this->mesh, mesh);
}

void InputOpenFoam::variables(Mesh &mesh)
{

}
