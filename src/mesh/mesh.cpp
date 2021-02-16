
#include "mesh.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include "basis/containers/serializededata.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/packing.h"
#include "wrappers/mpi/communication.h"

#include "input/builders/input.h"
#include "input/parsers/ansyscdb/ansyscdb.h"
#include "input/parsers/openfoam/openfoam.h"
#include "input/parsers/abaqus/abaqus.h"
#include "input/parsers/xdmf/xdmf.h"
#include "input/parsers/ensight/ensight.h"
#include "input/parsers/vtklegacy/vtklegacy.h"
#include "input/parsers/netgen/netgen.h"
#include "input/parsers/meshgenerator/meshgenerator.h"

#include "preprocessing/meshpreprocessing.h"
#include "store/statisticsstore.h"
#include "store/elementstore.h"
#include "store/nodestore.h"
#include "store/elementsregionstore.h"
#include "store/boundaryregionstore.h"
#include "store/contactinterfacestore.h"
#include "store/surfacestore.h"
#include "store/contactstore.h"
#include "store/fetidatastore.h"

#include "output/output.h"
#include "output/visualization/debug.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

using namespace espreso;

Element Mesh::edata[(int)Element::CODE::SIZE];

bool Mesh::convertDatabase()
{
	switch (info::ecf->input_type) {
	case ECF::INPUT_TYPE::EXTERNAL_FILE:
		return info::ecf->input.convert_database;
	default:
		return false;
	}
}

void Mesh::init()
{
	edata[static_cast<int>(Element::CODE::POINT1   )]   .init<Element::CODE::POINT1   >();
	edata[static_cast<int>(Element::CODE::LINE2    )]   .init<Element::CODE::LINE2    >();
	edata[static_cast<int>(Element::CODE::TRIANGLE3)]   .init<Element::CODE::TRIANGLE3>();
	edata[static_cast<int>(Element::CODE::SQUARE4  )]   .init<Element::CODE::SQUARE4  >();
	edata[static_cast<int>(Element::CODE::TETRA4   )]   .init<Element::CODE::TETRA4   >();
	edata[static_cast<int>(Element::CODE::PYRAMID5 )]   .init<Element::CODE::PYRAMID5 >();
	edata[static_cast<int>(Element::CODE::PRISMA6  )]   .init<Element::CODE::PRISMA6  >();
	edata[static_cast<int>(Element::CODE::HEXA8    )]   .init<Element::CODE::HEXA8    >();
	edata[static_cast<int>(Element::CODE::LINE3    )]   .init<Element::CODE::LINE3    >();
	edata[static_cast<int>(Element::CODE::TRIANGLE6)]   .init<Element::CODE::TRIANGLE6>();
	edata[static_cast<int>(Element::CODE::SQUARE8  )]   .init<Element::CODE::SQUARE8  >();
	edata[static_cast<int>(Element::CODE::TETRA10  )]   .init<Element::CODE::TETRA10  >();
	edata[static_cast<int>(Element::CODE::PYRAMID13)]   .init<Element::CODE::PYRAMID13>();
	edata[static_cast<int>(Element::CODE::PRISMA15 )]   .init<Element::CODE::PRISMA15 >();
	edata[static_cast<int>(Element::CODE::HEXA20   )]   .init<Element::CODE::HEXA20   >();

	info::mesh = new Mesh();
}

void Mesh::load()
{
	profiler::syncstart("load");
	MeshBuilder *data = NULL;
	switch (info::ecf->input_type) {
	case ECF::INPUT_TYPE::EXTERNAL_FILE:
		switch (info::ecf->input.format) {
		case InputConfiguration::FORMAT::ANSYS_CDB:      data = new AnsysCDBLoader     (info::ecf->input); break;
		case InputConfiguration::FORMAT::OPENFOAM:       data = new OpenFOAMLoader     (info::ecf->input); break;
		case InputConfiguration::FORMAT::ABAQUS:         data = new AbaqusLoader       (info::ecf->input); break;
		case InputConfiguration::FORMAT::XDMF:           data = new XDMFLoader         (info::ecf->input); break;
		case InputConfiguration::FORMAT::ENSIGHT:        data = new EnsightLoader      (info::ecf->input); break;
		case InputConfiguration::FORMAT::VTK_LEGACY:     data = new VTKLegacyLoader    (info::ecf->input); break;
		case InputConfiguration::FORMAT::NETGET: data = new NetgenNeutralLoader(info::ecf->input); break;
		}
		break;
	case ECF::INPUT_TYPE::GENERATOR:
	default:
		data = new MeshGenerator(info::ecf->generator);
	}

	data->load();
	data->build();

	delete data;

	profiler::syncend("load");
}

void Mesh::finish()
{
	delete info::mesh;
}

Mesh::Mesh()
: elements(new ElementStore()), nodes(new NodeStore()),
  FETIData(new FETIDataStore()),
  halo(new ElementStore()),
  surface(new SurfaceStore()), domainsSurface(new SurfaceStore()),
  contacts(new ContactStore()),

  output(new Output()),
  _withGUI(false)
{
	dimension = 0;
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D:
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
	case PhysicsConfiguration::TYPE::SHALLOW_WATER_2D:
		dimension = 2;
		break;
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D:
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		dimension = 3;
		break;
	}

	preferedDomains = info::ecf->input.decomposition.domains;
	if (preferedDomains == 0) {
		preferedDomains = info::env::OMP_NUM_THREADS; // TODO: set better value;
	}
	uniformDecomposition = false;
}

Mesh::~Mesh()
{
	// we need to delete output first in order to wait for unfinished output operations
	delete output;

	delete elements;
	delete nodes;
	delete FETIData;
	delete halo;
	delete surface;
	delete domainsSurface;
	delete contacts;

	for (size_t i = 0; i < contactInterfaces.size(); ++i) {
		delete contactInterfaces[i];
	}
	for (size_t i = 0; i < boundaryRegions.size(); ++i) {
		delete boundaryRegions[i];
	}
	for (size_t i = 0; i < elementsRegions.size(); ++i) {
		delete elementsRegions[i];
	}
}

ElementsRegionStore* Mesh::eregion(const std::string &name)
{
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (StringCompare::caseSensitiveEq(elementsRegions[r]->name, name)) {
			return elementsRegions[r];
		}
	}
	eslog::error("Unknown region of elements with name '%s'.\n", name.c_str());
	return NULL;
}

BoundaryRegionStore* Mesh::bregion(const std::string &name)
{
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (StringCompare::caseSensitiveEq(boundaryRegions[r]->name, name)) {
			return boundaryRegions[r];
		}
	}
	eslog::error("Unknown boundary region '%s'\n.", name.c_str());
	return NULL;
}

bool Mesh::onAllElements(const std::string &eregion) const
{
	return StringCompare::caseInsensitiveEq(eregion, "ALL_ELEMENTS");
}

bool Mesh::hasPhaseChange() const
{
	for (size_t m = 0; m < materials.size(); m++) {
		if (materials[m]->phase_change) {
			return true;
		}
	}
	return false;
}

void Mesh::setMaterials()
{
	materials.clear();
	std::map<std::string, int> matindex;
	for (auto mat = info::ecf->getPhysics()->materials.begin(); mat != info::ecf->getPhysics()->materials.end(); ++mat) {
		materials.push_back(&mat->second);
		matindex[mat->first] = materials.size() - 1;
	}

	for (auto mat = info::ecf->getPhysics()->material_set.begin(); mat != info::ecf->getPhysics()->material_set.end(); ++mat) {
		ElementsRegionStore *region = eregion(mat->first);
		if (matindex.find(mat->second) == matindex.end()) {
			eslog::globalerror("Unknown material '%s'.\n", mat->second.c_str());
		}
		int material = matindex.find(mat->second)->second;
		for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); ++e) {
			elements->material->datatarray()[*e] = material;
		}
	}
}

void Mesh::preprocess()
{
	profiler::syncstart("meshing");

	auto hasBEM = [] (const PhysicsConfiguration *physics) {
		for (auto it = physics->discretization.begin(); it != physics->discretization.end(); ++it) {
			if (it->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
				return true;
			}
		}
		return false;
	};

	auto forEachSteps = [&] (std::function<bool(const LoadStepSolverConfiguration &)> fnc) {
		bool ret = false;
		switch (info::ecf->physics) {
		case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D:
			for (auto step = info::ecf->thermo_elasticity_2d.load_steps_settings.begin(); step != info::ecf->thermo_elasticity_2d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second.heat_transfer);
				ret |= fnc(step->second.structural_mechanics);
			}
			break;
		case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D:
			for (auto step = info::ecf->thermo_elasticity_3d.load_steps_settings.begin(); step != info::ecf->thermo_elasticity_3d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second.heat_transfer);
				ret |= fnc(step->second.structural_mechanics);
			}
			break;
		case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
			for (auto step = info::ecf->heat_transfer_2d.load_steps_settings.begin(); step != info::ecf->heat_transfer_2d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
			for (auto step = info::ecf->heat_transfer_3d.load_steps_settings.begin(); step != info::ecf->heat_transfer_3d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
			for (auto step = info::ecf->structural_mechanics_2d.load_steps_settings.begin(); step != info::ecf->structural_mechanics_2d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
			for (auto step = info::ecf->structural_mechanics_3d.load_steps_settings.begin(); step != info::ecf->structural_mechanics_3d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		default:
			eslog::globalerror("Mesh: Not implemented physics.\n");
			exit(0);
		}
		return ret;
	};

	auto ntob = [&] (const std::string &rname, int dimension) {
		BoundaryRegionStore *region = bregion(rname);
		if (region->dimension == 0) {
			mesh::computeBoundaryElementsFromNodes(region, dimension);
		}
	};

	eslog::startln("MESH: PREPROCESSING STARTED", "MESHING");
	setMaterials();
	eslog::checkpointln("MESH: MATERIALS FILLED");

	if (info::ecf->input_type == ECF::INPUT_TYPE::EXTERNAL_FILE) {
		mesh::reclusterize();
		uniformDecomposition = false;
	}

	if (info::ecf->input_type == ECF::INPUT_TYPE::GENERATOR) {
		uniformDecomposition = info::ecf->generator.uniform_domains;
		if (!info::ecf->generator.uniform_clusters) {
			mesh::reclusterize();
			uniformDecomposition = false;
		}
	}
	profiler::synccheckpoint("reclusterize");

	esint msize;
	Communication::allReduce(&elements->size, &msize, 1, MPITools::getType<esint>().mpitype, MPI_MIN);
	if (msize == 0) {
		eslog::globalerror("ESPRESO quit computtion: process without any elements detected.\n");
	}

	if (hasBEM(info::ecf->getPhysics())) {
		// TODO: BEM does not always need separate regions
		info::ecf->input.decomposition.separate_materials = true;
		info::ecf->input.decomposition.separate_regions = true;
	}
	if (
			info::ecf->input.decomposition.separate_materials ||
			info::ecf->input.decomposition.separate_regions ||
			info::ecf->input.decomposition.separate_etypes) {
		uniformDecomposition = false;
	}
	mesh::partitiate(preferedDomains, uniformDecomposition);

	profiler::synccheckpoint("domain_decomposition");
	eslog::checkpointln("MESH: MESH DECOMPOSED");

	if (info::ecf->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D || info::ecf->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D) {
		const StructuralMechanicsConfiguration *sm;
		int dimension;
		if (info::ecf->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D) {
			sm = &info::ecf->structural_mechanics_2d;
			dimension = 1;
		}
		if (info::ecf->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D) {
			sm = &info::ecf->structural_mechanics_3d;
			dimension = 2;
		}

		for (auto ls = sm->load_steps_settings.begin(); ls != sm->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.normal_pressure.begin(); bc != ls->second.normal_pressure.end(); ++bc) {
				ntob(bc->first, dimension);
			}
		}
	}

	if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D || info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D) {
		const HeatTransferConfiguration *ht;
		int dimension;
		if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D) {
			ht = &info::ecf->heat_transfer_2d;
			dimension = 1;
		}
		if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D) {
			ht = &info::ecf->heat_transfer_3d;
			dimension = 2;
		}

		bool composeregion = false;
		for (auto ls = ht->load_steps_settings.begin(); ls != ht->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.heat_flow.begin(); bc != ls->second.heat_flow.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.heat_flux.begin(); bc != ls->second.heat_flux.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.convection.begin(); bc != ls->second.convection.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.diffuse_radiation.begin(); bc != ls->second.diffuse_radiation.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
		}

		if (composeregion) {
			eslog::checkpointln("MESH: BOUNDARY REGIONS COMPOSED");
		}
	}

	if (info::ecf->physics == PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D || info::ecf->physics == PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D) {
		const ThermoElasticityConfiguration *te;
		int dimension;
		if (info::ecf->physics == PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D) {
			te = &info::ecf->thermo_elasticity_2d;
			dimension = 1;
		}
		if (info::ecf->physics == PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D) {
			te = &info::ecf->thermo_elasticity_3d;
			dimension = 2;
		}

		bool composeregion = false;
		for (auto ls = te->load_steps_settings.begin(); ls != te->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.heat_transfer.heat_flow.begin(); bc != ls->second.heat_transfer.heat_flow.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.heat_transfer.heat_flux.begin(); bc != ls->second.heat_transfer.heat_flux.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.heat_transfer.convection.begin(); bc != ls->second.heat_transfer.convection.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.heat_transfer.diffuse_radiation.begin(); bc != ls->second.heat_transfer.diffuse_radiation.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.structural_mechanics.normal_pressure.begin(); bc != ls->second.structural_mechanics.normal_pressure.end(); ++bc) {
				composeregion = true;
				ntob(bc->first, dimension);
			}
		}

		if (composeregion) {
			eslog::checkpointln("MESH: BOUNDARY REGIONS COMPOSED");
		}
	}

	profiler::synccheckpoint("preprocess_boundary_regions");

	mesh::arrangeElementsRegions();
	profiler::synccheckpoint("arrange_element_regions");
	eslog::checkpointln("MESH: ELEMENT REGIONS ARRANGED");

	mesh::computeBodies();
	eslog::checkpointln("MESH: BODIES COMPUTED");

	if (info::ecf->input.contact_interfaces.size()) {
		mesh::computeBodiesSurface();
		mesh::computeWarpedNormals(surface);
		mesh::exchangeContactHalo();
		mesh::findCloseElements();
		mesh::computeContactInterface();
		mesh::arrangeContactInterfaces();
		profiler::synccheckpoint("compute_contact_interface");
		eslog::checkpointln("MESH: CONTACT INTERFACE COMPUTED");
	}

	mesh::arrangeBoundaryRegions();
	profiler::synccheckpoint("arrange_boundary_regions");
	eslog::checkpointln("MESH: BOUNDARY REGIONS ARRANGED");

	if (forEachSteps([] (const LoadStepSolverConfiguration &step) {
		return step.solver == LoadStepSolverConfiguration::SOLVER::FETI;
	})) {

		mesh::computeNodeDomainDistribution();
		mesh::computeLocalIndices();
		profiler::synccheckpoint("preprocess_domains");
		eslog::checkpointln("MESH: ELEMENTS DOMAIN INDICES COMPUTED");
	}

	if (info::ecf->getPhysics()->dimension == DIMENSION::D3 && (hasBEM(info::ecf->getPhysics()))) {
		mesh::computeDomainsSurface();
		mesh::triangularizeDomainSurface();
		profiler::synccheckpoint("preprocess_surface");
		eslog::checkpointln("MESH: DOMAIN SURFACE COMPUTED");
	}

	if (info::ecf->getPhysics()->dimension == DIMENSION::D3 && _withGUI) {
		mesh::computeRegionsSurface();
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			mesh::triangularizeSurface(elementsRegions[r]->surface);
		}
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			mesh::triangularizeBoundary(boundaryRegions[r]);
		}
		eslog::checkpointln("MESH: REGION SURFACE COMPUTED");
	}

	if (info::ecf->getPhysics()->dimension == DIMENSION::D3 && info::ecf->output.format == OutputConfiguration::FORMAT::STL_SURFACE) {
		mesh::computeBodiesSurface();
		mesh::triangularizeSurface(surface);
		eslog::checkpointln("MESH: BODIES SURFACE COMPUTED");
	}

	if (forEachSteps([] (const LoadStepSolverConfiguration &step) {
		return step.solver == LoadStepSolverConfiguration::SOLVER::FETI && step.feti.method == FETIConfiguration::METHOD::HYBRID_FETI && step.feti.B0_type == FETIConfiguration::B0_TYPE::KERNELS;
	})) {

		mesh::computeDomainDual();
		profiler::synccheckpoint("compute_domain_dual");
	}

	if (info::ecf->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D) {
		if (info::ecf->structural_mechanics_3d.discretization.size() && info::ecf->structural_mechanics_3d.discretization.begin()->second == PhysicsConfiguration::DISCRETIZATION::FEM_TDNNS) {
			mesh::computeElementsEdgeNeighbors();
			profiler::synccheckpoint("compute_edge_neighbors");
			eslog::checkpointln("MESH: EDGE NEIGHBORS COMPUTED");
		}
	}

	DebugOutput::mesh();
	profiler::syncend("meshing");
	eslog::endln("MESH: PREPROCESSING FINISHED");
}

void Mesh::duplicate()
{
	eslog::startln("MESH: CREATE DUPLICATED INSTANCES", "DUPLICATION");

	size_t packedSize = 0;

	if (info::mpi::irank == 0) {
		packedSize += utils::packedSize(dimension);
		packedSize += utils::packedSize(preferedDomains);
		packedSize += utils::packedSize(uniformDecomposition);

		packedSize += elements->packedFullSize();
		packedSize += nodes->packedFullSize();

		packedSize += utils::packedSize(elementsRegions.size());
		for (size_t i = 0; i < elementsRegions.size(); i++) {
			packedSize += elementsRegions[i]->packedFullSize();
		}
		packedSize += utils::packedSize(boundaryRegions.size());
		for (size_t i = 0; i < boundaryRegions.size(); i++) {
			packedSize += boundaryRegions[i]->packedFullSize();
		}
		packedSize += utils::packedSize(contactInterfaces.size());
		for (size_t i = 0; i < contactInterfaces.size(); i++) {
			packedSize += contactInterfaces[i]->packedFullSize();
		}

		packedSize += FETIData->packedFullSize();
		packedSize += halo->packedFullSize();

		packedSize += surface->packedFullSize();
		packedSize += domainsSurface->packedFullSize();
		packedSize += contacts->packedFullSize();

		packedSize += utils::packedSize(neighbors);
		packedSize += utils::packedSize(neighborsWithMe);
		packedSize += utils::packedSize(_withGUI);
	}

	Communication::broadcast(&packedSize, sizeof(size_t), MPI_BYTE, 0, MPITools::instances);
	char *buffer = new char[packedSize];

	if (info::mpi::irank == 0) {
		char *p = buffer;
		utils::pack(dimension, p);
		utils::pack(preferedDomains, p);
		utils::pack(uniformDecomposition, p);

		elements->packFull(p);
		nodes->packFull(p);

		utils::pack(elementsRegions.size(), p);
		for (size_t i = 0; i < elementsRegions.size(); i++) {
			elementsRegions[i]->packFull(p);
		}
		utils::pack(boundaryRegions.size(), p);
		for (size_t i = 0; i < boundaryRegions.size(); i++) {
			boundaryRegions[i]->packFull(p);
		}
		utils::pack(contactInterfaces.size(), p);
		for (size_t i = 0; i < contactInterfaces.size(); i++) {
			contactInterfaces[i]->packFull(p);
		}

		FETIData->packFull(p);
		halo->packFull(p);

		surface->packFull(p);
		domainsSurface->packFull(p);
		contacts->packFull(p);

		utils::pack(neighbors, p);
		utils::pack(neighborsWithMe, p);
		utils::pack(_withGUI, p);
	}

	eslog::checkpoint("MESH: MESH PACKED");
	eslog::param("size[MB]", packedSize);
	eslog::ln();

	Communication::broadcast(buffer, packedSize, MPI_CHAR, 0, MPITools::instances);

	eslog::checkpointln("MESH: PACKED DATA BROADCASTED");

	if (info::mpi::irank != 0) {
		for (size_t i = 0; i < elementsRegions.size(); i++) {
			delete elementsRegions[i];
		}
		elementsRegions.clear();
		for (size_t i = 0; i < boundaryRegions.size(); i++) {
			delete boundaryRegions[i];
		}
		boundaryRegions.clear();

		const char *p = buffer;
		utils::unpack(dimension, p);
		utils::unpack(preferedDomains, p);
		utils::unpack(uniformDecomposition, p);

		elements->unpackFull(p);
		nodes->unpackFull(p);

		size_t size;
		utils::unpack(size, p);
		for (size_t i = 0; i < size; i++) {
			elementsRegions.push_back(new ElementsRegionStore(p));
		}
		utils::unpack(size, p);
		for (size_t i = 0; i < size; i++) {
			boundaryRegions.push_back(new BoundaryRegionStore(p));
		}
		utils::unpack(size, p);
		for (size_t i = 0; i < size; i++) {
			contactInterfaces.push_back(new ContactInterfaceStore(p));
		}

		FETIData->unpackFull(p);
		halo->unpackFull(p);

		surface->unpackFull(p);
		domainsSurface->unpackFull(p);
		contacts->unpackFull(p);

		utils::unpack(neighbors, p);
		utils::unpack(neighborsWithMe, p);
		utils::unpack(_withGUI, p);

		setMaterials();
	}

	delete[] buffer;

	eslog::endln("MESH: DUPLICATION FINISHED");
}

void Mesh::toBuffer()
{
	for (size_t i = 0; i < elements->data.size(); ++i) {
		if (elements->data[i]->name.size()) {
			elements->data[i]->toBuffer();
		}
	}
	for (size_t i = 0; i < nodes->data.size(); ++i) {
		if (nodes->data[i]->name.size()) {
			nodes->data[i]->toBuffer();
		}
	}
}

void Mesh::printMeshStatistics()
{
	size_t namesize = 56;

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (namesize < boundaryRegions[r]->name.size() + 16) {
			namesize = boundaryRegions[r]->name.size() + 16;
		}
	}

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (namesize < elementsRegions[r]->name.size() + 16) {
			namesize = elementsRegions[r]->name.size() + 16;
		}
	}

	auto ename = [] (int code) -> std::string {
		switch (static_cast<Element::CODE>(code)) {

		case Element::CODE::POINT1: return "POINT1";

		case Element::CODE::LINE2: return "LINE2";

		case Element::CODE::TRIANGLE3: return "TRIANGLE3";
		case Element::CODE::SQUARE4: return "SQUARE4";

		case Element::CODE::TETRA4: return "TETRA4";
		case Element::CODE::PYRAMID5: return "PYRAMID5";
		case Element::CODE::PRISMA6: return "PRISMA6";
		case Element::CODE::HEXA8: return "HEXA8";

		case Element::CODE::LINE3: return "LINE3";

		case Element::CODE::TRIANGLE6: return "TRIANGLE6";
		case Element::CODE::SQUARE8: return "SQUARE8";

		case Element::CODE::TETRA10: return "TETRA10";
		case Element::CODE::PYRAMID13: return "PYRAMID13";
		case Element::CODE::PRISMA15: return "PRISMA15";
		case Element::CODE::HEXA20: return "HEXA20";

		default:
			eslog::internalFailure("unknown element code.\n");
			return "";
		}
	};

	size_t nelements = 0, nelementsnodes = 0;
	size_t fregs = 0, nfaces = 0, nfacenodes = 0;
	size_t eregs = 0, nedges = 0, nedgenodes = 0;
	size_t nregs = 0, nnodes = 0;

	for (size_t r = 1; r < elementsRegions.size(); r++) {
		nelements += elementsRegions[r]->totalsize; nelementsnodes += elementsRegions[r]->nodeInfo.totalSize;
	}
	for (size_t r = 1; r < boundaryRegions.size(); r++) {
		switch (boundaryRegions[r]->dimension) {
		case 2: ++fregs; nfaces += boundaryRegions[r]->totalsize; nfacenodes += boundaryRegions[r]->nodeInfo.totalSize; break;
		case 1: ++eregs; nedges += boundaryRegions[r]->totalsize; nedgenodes += boundaryRegions[r]->nodeInfo.totalSize; break;
		case 0: ++nregs; nnodes += boundaryRegions[r]->nodeInfo.totalSize; break;
		default: break;
		}
	}

	esint ecountTotal = 0, scountTotal = 0;
	for (size_t b = 0; b < elementsRegions[0]->bodies.size(); ++b) {
		ecountTotal += elementsRegions[0]->bodyElements[b];
		scountTotal += elementsRegions[0]->bodyFaces[b];
	}

	switch (info::ecf->output.logger) {
	case OutputConfiguration::LOGGER::USER:
		eslog::info(" ====================================== MESH STATISTICS ====================================== \n");
		eslog::info(" ============================================================================================= \n");
		eslog::info("  %s%*s : %16s %16s %16s\n", "GEOMETRY", namesize - 26, " ", "MIN", "MAX", "LENGTH");
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		eslog::info("  %*s : %16f %16f %16f\n", namesize - 18, "X-COORDINATES", elementsRegions.front()->nodeInfo.min.x, elementsRegions.front()->nodeInfo.max.x, elementsRegions.front()->nodeInfo.max.x - elementsRegions.front()->nodeInfo.min.x);
		eslog::info("  %*s : %16f %16f %16f\n", namesize - 18, "Y-COORDINATES", elementsRegions.front()->nodeInfo.min.y, elementsRegions.front()->nodeInfo.max.y, elementsRegions.front()->nodeInfo.max.y - elementsRegions.front()->nodeInfo.min.y);
		eslog::info("  %*s : %16f %16f %16f\n", namesize - 18, "Z-COORDINATES", elementsRegions.front()->nodeInfo.min.z, elementsRegions.front()->nodeInfo.max.z, elementsRegions.front()->nodeInfo.max.z - elementsRegions.front()->nodeInfo.min.z);
		eslog::info(" ============================================================================================= \n");

		eslog::info(" %*s : %16s %16s\n", namesize, "REGION NAME", "ELEMENTS", "NODES");
		eslog::info(" ============================================================================================= \n");

		eslog::info("  %s%*s : %16s %16s\n", "TOTAL NUMBER OF NODES", namesize - 22, " ", " ", Parser::stringwithcommas(nodes->uniqInfo.totalSize).c_str());
		eslog::info("  %s%*s : %16s\n", "TOTAL NUMBER OF ELEMENTS", namesize - 25, " ", Parser::stringwithcommas(elements->totalSize).c_str());

		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (elements->ecounters[etype]) {
				eslog::info(" %*s : %16s\n", namesize, ename(etype).c_str(), Parser::stringwithcommas(elements->ecounters[etype]).c_str());
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n\n");

		eslog::info("  %s [%3ld]%*s : %16s %16s\n", "ELEMS REGIONS SIZES", elementsRegions.size() - 1, namesize - 26, " ", Parser::stringwithcommas(nelements).c_str(), Parser::stringwithcommas(nelementsnodes).c_str());
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(elementsRegions[r]->name, "NAMELESS_ELEMENT_SET")) {
				eslog::warning(" %*s : %16s %16s\n", namesize, elementsRegions[r]->name.c_str(), Parser::stringwithcommas(elementsRegions[r]->totalsize).c_str(), Parser::stringwithcommas(elementsRegions[r]->nodeInfo.totalSize).c_str());
			} else {
				eslog::info(" %*s : %16s %16s\n", namesize, elementsRegions[r]->name.c_str(), Parser::stringwithcommas(elementsRegions[r]->totalsize).c_str(), Parser::stringwithcommas(elementsRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "FACES REGIONS SIZES", fregs, namesize - 26, " ", Parser::stringwithcommas(nfaces).c_str(), Parser::stringwithcommas(nfacenodes).c_str());
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->dimension == 2) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), Parser::stringwithcommas(boundaryRegions[r]->totalsize).c_str(), Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "EDGES REGIONS SIZES", eregs, namesize - 26, " ", Parser::stringwithcommas(nedges).c_str(), Parser::stringwithcommas(nedgenodes).c_str());
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->dimension == 1) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), Parser::stringwithcommas(boundaryRegions[r]->totalsize).c_str(), Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "NODES REGIONS SIZES", nregs, namesize - 26, " ", " ", Parser::stringwithcommas(nnodes).c_str());
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->dimension == 0) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), " ", Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info(" ============================================================================================= \n");

		eslog::info("  BODY STATISTICS %22s : %16s %16s %16s\n", "ELEMENTS REGION", "ELEMENTS", "FACES", "PROPORTION");
		eslog::info(" ============================================================================================= \n");
		eslog::info("  %38s : %16s %16s %9d BODIES\n", elementsRegions[0]->name.c_str(), Parser::stringwithcommas(ecountTotal).c_str(), Parser::stringwithcommas(scountTotal).c_str(), elements->bodiesTotalSize);
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		for (size_t r = 1; r < elementsRegions.size(); r++) {
			for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
				if (elementsRegions[r]->bodyElements[b] == elementsRegions[0]->bodyElements[elementsRegions[r]->bodies[b]] && elementsRegions[r]->bodies.size() == 1) {
					eslog::info("  %38s : %16s %16s %16s\n", b ? "" : elementsRegions[r]->name.c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyElements[b]).c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyFaces[b]).c_str(), "BODY=REGION");
				} else {
					double ratio = (double)elementsRegions[r]->bodyElements[b] / elementsRegions[0]->bodyElements[elementsRegions[r]->bodies[b]];
					eslog::info("  %38s : %16s %16s %16f\n", b ? "" : elementsRegions[r]->name.c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyElements[b]).c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyFaces[b]).c_str(), ratio);
				}
			}
		}
		eslog::info(" ============================================================================================= \n");

		if (info::ecf->input.contact_interfaces.size()) {
			eslog::info("  CONTACT INTERFACES %16s %12s : %27s %12s\n", "FACES", "AREA", "FACES", "AREA");
			eslog::info(" ============================================================================================= \n");
		}
		for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
			std::string interface = " - - - - " + it->first;
			switch (it->second.detection) {
			case ContactInterfaceConfiguration::DETECTION::ALL_BODIES: interface += ", ALL_BODIES"; break;
			case ContactInterfaceConfiguration::DETECTION::BODY_LIST: interface += ", BODY_LIST"; break;
			case ContactInterfaceConfiguration::DETECTION::CONTACT_PAIR: interface += ", CONTACT_PAIR"; break;
			}
			switch (it->second.criterion) {
			case ContactInterfaceConfiguration::CRITERION::BOUND: interface += ", BOUND"; break;
			case ContactInterfaceConfiguration::CRITERION::GAP: interface += ", GAP=" + std::to_string(it->second.gap); break;
			}
			interface += " - - - - ";
			eslog::info(" %*s %*s \n", 36, " ", namesize, interface.c_str());
			for (size_t i = 0; i < it->second.found_interfaces.size(); ++i) {
				auto &iface = contacts->interfaces[contactInterfaces[it->second.found_interfaces[i]]->interfaceIndex];
				std::vector<std::string> names = Parser::split(contactInterfaces[it->second.found_interfaces[i]]->name, "-");
				if (names.size() <= 2) {
					names.push_back("NAMELESS_BODY");
				}
				if (names.size() <= 3) {
					names.push_back("NAMELESS_BODY");
				}
				eslog::info("                                                                                               \n");
				eslog::info(" %*s %30s : %40s \n", 18, "", names[2].c_str(), names[3].c_str());
				eslog::info("          %27s %12f : %27s %12f\n", Parser::stringwithcommas(iface.from.faces).c_str(), iface.from.area, Parser::stringwithcommas(iface.to.faces).c_str(), iface.to.area);
			}
		}
		eslog::info(" ============================================================================================= \n");
		break;
	case OutputConfiguration::LOGGER::PARSER:
		eslog::info(" ====================================== MESH STATISTICS ====================================== \n");
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			eslog::info("mesh: region=%s, dimension=%d, elements=%d, nodes=%d\n", elementsRegions[r]->name.c_str(), dimension, elementsRegions[r]->totalsize, elementsRegions[r]->nodeInfo.totalSize);
		}
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			eslog::info("mesh: region=%s, dimension=%d, elements=%d, nodes=%d\n", boundaryRegions[r]->name.c_str(), boundaryRegions[r]->dimension, boundaryRegions[r]->totalsize, boundaryRegions[r]->nodeInfo.totalSize);
		}
		eslog::info("mesh: region=ALL_ELEMENTS, bodies=%d\n", elements->bodiesTotalSize);
		for (size_t r = 1; r < elementsRegions.size(); r++) {
			for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
				eslog::info("mesh: region=%s, b-elements=%d, b-faces=%d\n", elementsRegions[r]->name.c_str(), elementsRegions[r]->bodyElements[b], elementsRegions[r]->bodyFaces[b]);
			}
		}
		for (auto it = contactInterfaces.begin(); it != contactInterfaces.end(); ++it) {
			eslog::info("mesh: region=%s, s-faces=%d, s-area=%.5f, d-faces=%d, d-area=%.5f\n", (*it)->name.c_str(),
					contacts->interfaces[(*it)->interfaceIndex].from.faces, contacts->interfaces[(*it)->interfaceIndex].from.area,
					contacts->interfaces[(*it)->interfaceIndex].to.faces, contacts->interfaces[(*it)->interfaceIndex].to.area);
		}
		break;
	}
}

void Mesh::printDecompositionStatistics()
{
	esint cluster     = 0;
	esint domain      = 1;
	esint elements    = 2;
	esint nodes       = 3;
	esint neighbors   = 4;
	esint elPerDomain = 5;
//	esint nPerDomain  = 6;
	esint total       = 7;
	std::vector<esint> min(total, std::numeric_limits<esint>::max()), max(total), sum(total), gmin(total, 1), gmax(total, 1), gsum(total, 1);

	min[cluster] = max[cluster] = sum[cluster] = this->elements->nclusters;
	min[domain] = max[domain] = sum[domain] = this->elements->ndomains;
	min[elements] = max[elements] = sum[elements] = this->elements->size;
	min[nodes] = max[nodes] = sum[nodes] = this->nodes->size;
	min[neighbors] = max[neighbors] = sum[neighbors] = this->neighbors.size();

	for (esint d = 0; d < this->elements->ndomains; ++d) {
		min[elPerDomain] = std::min(min[elPerDomain], this->elements->elementsDistribution[d + 1] - this->elements->elementsDistribution[d]);
		max[elPerDomain] = std::max(max[elPerDomain], this->elements->elementsDistribution[d + 1] - this->elements->elementsDistribution[d]);
		sum[elPerDomain] += this->elements->elementsDistribution[d + 1] - this->elements->elementsDistribution[d];

//		esint nodes = this->nodes->uniqInfo.totalSize;
//		eslog::warning("DOMAIN INTERVALS DEPENDENT CODE\n");
//		esint nodes = 0;
//		for (size_t i = 0; i < this->nodes->dintervals[d].size(); ++i) {
//			nodes += this->nodes->dintervals[d][i].end - this->nodes->dintervals[d][i].begin;
//		}
//		min[nPerDomain] = std::min(min[nPerDomain], nodes);
//		max[nPerDomain] = std::max(max[nPerDomain], nodes);
//		sum[nPerDomain] += nodes;
	}

	Communication::reduce(min.data(), gmin.data(), total, MPITools::getType<esint>().mpitype, MPI_MIN, 0);
	Communication::reduce(max.data(), gmax.data(), total, MPITools::getType<esint>().mpitype, MPI_MAX, 0);
	Communication::reduce(sum.data(), gsum.data(), total, MPITools::getType<esint>().mpitype, MPI_SUM, 0);

	double imbalance;
	eslog::info("                                                                                               \n");
	eslog::info(" ================================== DECOMPOSITION STATISTICS ================================= \n");
	eslog::info("                              MIN                 MAX                 SUM           IMBALANCE \n");
	eslog::info("  CLUSTERS  : %19s %19s %19s", Parser::stringwithcommas(gmin[cluster]).c_str(), Parser::stringwithcommas(gmax[cluster]).c_str(), Parser::stringwithcommas(gsum[cluster]).c_str());
	if ((imbalance = (double)gmax[cluster] / gmin[cluster]) > 2) {
		eslog::warning(" %19.3f\n", imbalance);
	} else {
		eslog::info(" %19.3f\n", imbalance);
	}
	eslog::info("  DOMAINS   : %19s %19s %19s", Parser::stringwithcommas(gmin[domain]).c_str(), Parser::stringwithcommas(gmax[domain]).c_str(), Parser::stringwithcommas(gsum[domain]).c_str());
	if ((imbalance = (double)gmax[domain] / gmin[domain]) > 2) {
		eslog::warning(" %19.3f\n", imbalance);
	} else {
		eslog::info(" %19.3f\n", imbalance);
	}
	eslog::info("  ELEMENTS  : %19s %19s %19s", Parser::stringwithcommas(gmin[elements]).c_str(), Parser::stringwithcommas(gmax[elements]).c_str(), Parser::stringwithcommas(gsum[elements]).c_str());
	if ((imbalance = (double)gmax[elements] / gmin[elements]) > 2) {
		eslog::warning(" %19.3f\n", imbalance);
	} else {
		eslog::info(" %19.3f\n", imbalance);
	}
	eslog::info("  NODES     : %19s %19s %19s", Parser::stringwithcommas(gmin[nodes]).c_str(), Parser::stringwithcommas(gmax[nodes]).c_str(), Parser::stringwithcommas(gsum[nodes]).c_str());
	if ((imbalance = (double)gmax[nodes] / gmin[nodes]) > 2) {
		eslog::warning(" %19.3f\n", imbalance);
	} else {
		eslog::info(" %19.3f\n", imbalance);
	}
	eslog::info("  NEIGHBORS : %19s %19s %19s", Parser::stringwithcommas(gmin[neighbors]).c_str(), Parser::stringwithcommas(gmax[neighbors]).c_str(), Parser::stringwithcommas(gsum[neighbors]).c_str());
	eslog::info(" %19.3f\n", gmin[neighbors] ? (double)gmax[neighbors] / gmin[neighbors] : 0);
	eslog::info("  ELEMS/DOM : %19s %19s %19s", Parser::stringwithcommas(gmin[elPerDomain]).c_str(), Parser::stringwithcommas(gmax[elPerDomain]).c_str(), Parser::stringwithcommas(gsum[elPerDomain]).c_str());
	if ((imbalance = (double)gmax[elPerDomain] / gmin[elPerDomain]) > 2) {
		eslog::warning(" %19.3f\n", imbalance);
	} else {
		eslog::info(" %19.3f\n", imbalance);
	}
//	eslog::info("  NODES/DOM : %13d %13d %13d", gmin[nPerDomain], gmax[nPerDomain], gsum[nPerDomain]);
//	if ((imbalance = (double)gmax[nPerDomain] / gmin[nPerDomain]) > 2) {
//		eslog::warning(" %13.3f\n", imbalance);
//	} else {
//		eslog::info(" %13.3f\n", imbalance);
//	}
	eslog::info(" ============================================================================================= \n");
}

