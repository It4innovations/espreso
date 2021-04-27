
#include "mesio.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/utilities/sysutils.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"

#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/bodystore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "wrappers/mpi/communication.h"

#include <cstdio>

struct MESIOData {
	MESIOData()
	{
		etypes.reserve(espreso::info::mesh->elements->process.size);
		for (size_t e = 0; e < espreso::info::mesh->elements->epointers->datatarray().size(); ++e) {
			etypes.push_back(static_cast<int>(espreso::info::mesh->elements->epointers->datatarray()[e]->code));
		}
		domains.reserve(espreso::info::mesh->elements->process.size);
		for (size_t i = 1; i < espreso::info::mesh->domains->elements.size(); ++i) {
			domains.resize(espreso::info::mesh->domains->elements[i], espreso::info::mesh->domains->offset + i - 1);
		}
		rtypes.resize(espreso::info::mesh->boundaryRegions.size());
		for (size_t i = 0; i < espreso::info::mesh->boundaryRegions.size(); ++i) {
			if (espreso::info::mesh->boundaryRegions[i]->dimension) {
				rtypes[i].reserve(espreso::info::mesh->boundaryRegions[i]->epointers->datatarray().size());
				for (size_t e = 0; e < espreso::info::mesh->boundaryRegions[i]->epointers->datatarray().size(); ++e) {
					rtypes[i].push_back(static_cast<int>(espreso::info::mesh->boundaryRegions[i]->epointers->datatarray()[e]->code));
				}
			}
		}
	}

	std::vector<esint> etypes, domains;
	std::vector<std::vector<esint> > rtypes;
};

using namespace espreso;

void MESIOInit(
	MPI_Comm		comm,
	int				verbosity)
{
	info::env::set();
	info::mpi::init(comm);
	MPITools::init();
	eslog::init(new Logger<ProgressTerminalLogger>);
	if (utils::exists("espreso.ecf")) {
		ECF::init("espreso.ecf");
	} else {
		ECF::init();
	}
	for (int i = 0; i < eslog::logger->size; ++i) {
		eslog::logger->args[i]->verbosity = verbosity;
	}
	Mesh::init();

	eslog::startln("MESIO: INITIALIZED", "MESIO");
}

void MESIOFinalize()
{
	Mesh::finish();
	ECF::finish();
	MPITools::finish();
	eslog::endln("MESIO: FINISHED");
}

void MESIOLoad(
	MESIO*			mesio,
	MESIOFormat		format,
	const char*		path,
	MESIODecomposer	decomposer,
	int				domains)
{
	switch (format) {
	case MESIO_ANSYS: info::ecf->input.format = InputConfiguration::FORMAT::ANSYS_CDB; break;
	case MESIO_ENSIGHT: info::ecf->input.format = InputConfiguration::FORMAT::ENSIGHT; break;
	case MESIO_VTK_LEGACY: info::ecf->input.format = InputConfiguration::FORMAT::VTK_LEGACY; break;
	case MESIO_XDMF: info::ecf->input.format = InputConfiguration::FORMAT::XDMF; break;
	}
	info::ecf->input.path = path;
	switch (decomposer) {
	case MESIO_NONE: info::ecf->input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::NONE; break;
	case MESIO_METIS: info::ecf->input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::METIS; break;
	case MESIO_PARMETIS: info::ecf->input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::PARMETIS; break;
	case MESIO_PTSCOTCH: info::ecf->input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::PTSCOTCH; break;
	case MESIO_HILBERT_CURVE: info::ecf->input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE; break;
	}
	info::ecf->input.decomposition.domains = domains;
	info::mesh->preferedDomains = domains;

	Mesh::load();
	info::mesh->preprocess();

	*mesio = new MESIOData();
}

void MESIONodes(
	MESIO           mesio,
	MESIOInt*       nhalo,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      ids,
	MESIOInt**      position,
	MESIOReal**     coordinates)
{
	*nhalo = info::mesh->nodes->uniqInfo.nhalo;
	*offset = info::mesh->nodes->uniqInfo.offset;
	*size = info::mesh->nodes->uniqInfo.size;
	*totalSize = info::mesh->nodes->uniqInfo.totalSize;

	*ids = info::mesh->nodes->IDs->datatarray().data();
	*position = info::mesh->nodes->uniqInfo.position.data();
	*coordinates = static_cast<double*>(&info::mesh->nodes->coordinates->datatarray()[0].x);
}

void MESIONodesRanks(
	MESIO           mesio,
	MESIOInt**      rankDistribution,
	int**           rankData)
{
	*rankDistribution = info::mesh->nodes->ranks->boundarytarray().data();
	*rankData = info::mesh->nodes->ranks->datatarray().data();
}

void MESIONodesDomains(
	MESIO           mesio,
	MESIOInt**      domainDistribution,
	int**           domainData)
{
	*domainDistribution = info::mesh->nodes->domains->boundarytarray().data();
	*domainData = info::mesh->nodes->domains->datatarray().data();
}

void MESIOElements(
	MESIO           mesio,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      type,
	MESIOInt**      enodesDistribution,
	MESIOInt**      enodesData)
{
	*offset = info::mesh->elements->process.offset;
	*size = info::mesh->elements->process.size;
	*totalSize = info::mesh->elements->process.totalSize;
	*type = mesio->etypes.data();
	*enodesDistribution = info::mesh->elements->nodes->boundarytarray().data();
	*enodesData = info::mesh->elements->nodes->datatarray().data();
}

void MESIOElementsDomains(
	MESIO           mesio,
	MESIOInt**      domains)
{
	*domains = mesio->domains.data();
}

void MESIOElementsMaterials(
	MESIO           mesio,
	int**           material)
{
	*material = info::mesh->elements->material->datatarray().data();
}

void MESIOElementsBodies(
	MESIO           mesio,
	int*            bodies,
	int**           body)
{
	*bodies = info::mesh->bodies->totalSize;
	*body = info::mesh->elements->body->datatarray().data();
}

void MESIOElementsNeighbors(
	MESIO           mesio,
	MESIOInt**      neighborDistribution,
	MESIOInt**      neighborData)
{
	*neighborDistribution = info::mesh->elements->faceNeighbors->boundarytarray().data();
	*neighborData = info::mesh->elements->faceNeighbors->datatarray().data();
}

void MESIOElementsCounters(
	MESIO           mesio,
	MESIOInt        etype,
	MESIOInt*       offset,
	MESIOInt*       totalSize)
{
	*offset = info::mesh->elementsRegions.front()->processPerCode[etype].offset;
	*totalSize = info::mesh->elementsRegions.front()->processPerCode[etype].totalSize;
}

int MESIOElementsRegions(
	MESIO           mesio)
{
	return info::mesh->elementsRegions.size();
}

void MESIOElementsRegion(
	MESIO           mesio,
	MESIOInt        region,
	const char**    name,
	MESIOInt*       size,
	MESIOInt**      elements)
{
	*name = info::mesh->elementsRegions[region]->name.c_str();
	*size = info::mesh->elementsRegions[region]->elements->datatarray().size();
	*elements = info::mesh->elementsRegions[region]->elements->datatarray().data();
}

void MESIOElementsRegionNodes(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt*       nhalo,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      nodes,
	MESIOInt**      position)
{
	*nhalo = info::mesh->elementsRegions[region]->nodeInfo.nhalo;
	*offset = info::mesh->elementsRegions[region]->nodeInfo.offset;
	*size = info::mesh->elementsRegions[region]->nodeInfo.size;
	*totalSize = info::mesh->elementsRegions[region]->nodeInfo.totalSize;
	*nodes = info::mesh->elementsRegions[region]->nodes->datatarray().data();
	*position = info::mesh->elementsRegions[region]->nodeInfo.position.data();
}

void MESIOElementsRegionCounters(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt        etype,
	MESIOInt*       offset,
	MESIOInt*       totalSize)
{
	*offset = info::mesh->elementsRegions[region]->processPerCode[etype].offset;
	*totalSize = info::mesh->elementsRegions[region]->processPerCode[etype].totalSize;
}

int MESIOBoundaryRegions(
	MESIO           mesio)
{
	return info::mesh->boundaryRegions.size();
}

void MESIOBoundaryRegion(
	MESIO           mesio,
	MESIOInt        region,
	const char**    name,
	MESIOInt*       dimension,
	MESIOInt*       size,
	MESIOInt**      type,
	MESIOInt**      parent,
	MESIOInt**      elementDistribution,
	MESIOInt**      elementData)
{
	*name = info::mesh->boundaryRegions[region]->name.c_str();
	*dimension = info::mesh->boundaryRegions[region]->dimension;
	if (*dimension) {
		*size = info::mesh->boundaryRegions[region]->process.size;
		*type = mesio->rtypes[region].data();
		*parent = info::mesh->boundaryRegions[region]->emembership->datatarray().data();
		*elementDistribution = info::mesh->boundaryRegions[region]->procNodes->boundarytarray().data();
		*elementData = info::mesh->boundaryRegions[region]->procNodes->datatarray().data();
	} else {
		*size = 0;
		*type = *parent = *elementDistribution = *elementData = NULL;
	}
}

void MESIOBoundaryRegionNodes(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt*       nhalo,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      nodes,
	MESIOInt**      position)
{
	*nhalo = info::mesh->boundaryRegions[region]->nodeInfo.nhalo;
	*offset = info::mesh->boundaryRegions[region]->nodeInfo.offset;
	*size = info::mesh->boundaryRegions[region]->nodeInfo.size;
	*totalSize = info::mesh->boundaryRegions[region]->nodeInfo.totalSize;
	*nodes = info::mesh->boundaryRegions[region]->nodes->datatarray().data();
	*position = info::mesh->boundaryRegions[region]->nodeInfo.position.data();
}

void MESIOBoundaryRegionCounters(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt        etype,
	MESIOInt*       offset,
	MESIOInt*       totalSize)
{
	*offset = info::mesh->boundaryRegions[region]->processPerCode[etype].offset;
	*totalSize = info::mesh->boundaryRegions[region]->processPerCode[etype].totalSize;
}
