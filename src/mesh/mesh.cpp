
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
#include "input/parsers/gmsh/gmsh.h"
#include "input/parsers/nglib/nglib.h"
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
		case InputConfiguration::FORMAT::NETGET:         data = new NetgenNeutralLoader(info::ecf->input); break;
		case InputConfiguration::FORMAT::GMSH:           data = new GMSHGenerator      (info::ecf->input); break;
		case InputConfiguration::FORMAT::NGLIB:          data = new NGLibGenerator     (info::ecf->input); break;
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
  _omitClusterization(false),
  _omitDecomposition(false),
  _withGUI(false),
  _withFETI(false),
  _withBEM(false),
  _withEdgeDual(false)
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

void Mesh::analyze()
{
	// check FETI
	switch (info::ecf->physics) {
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_2D:
		for (auto step = info::ecf->thermo_elasticity_2d.load_steps_settings.begin(); step != info::ecf->thermo_elasticity_2d.load_steps_settings.end(); ++step) {
			_withFETI |= step->second.heat_transfer.solver == LoadStepSolverConfiguration::SOLVER::FETI;
			_withFETI |= step->second.structural_mechanics.solver == LoadStepSolverConfiguration::SOLVER::FETI;
		}
		break;
	case PhysicsConfiguration::TYPE::THERMO_ELASTICITY_3D:
		for (auto step = info::ecf->thermo_elasticity_3d.load_steps_settings.begin(); step != info::ecf->thermo_elasticity_3d.load_steps_settings.end(); ++step) {
			_withFETI |= step->second.heat_transfer.solver == LoadStepSolverConfiguration::SOLVER::FETI;
			_withFETI |= step->second.structural_mechanics.solver == LoadStepSolverConfiguration::SOLVER::FETI;
		}
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		for (auto step = info::ecf->heat_transfer_2d.load_steps_settings.begin(); step != info::ecf->heat_transfer_2d.load_steps_settings.end(); ++step) {
			_withFETI |= step->second.solver == LoadStepSolverConfiguration::SOLVER::FETI;
		}
		break;
	case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		for (auto step = info::ecf->heat_transfer_3d.load_steps_settings.begin(); step != info::ecf->heat_transfer_3d.load_steps_settings.end(); ++step) {
			_withFETI |= step->second.solver == LoadStepSolverConfiguration::SOLVER::FETI;
		}
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
		for (auto step = info::ecf->structural_mechanics_2d.load_steps_settings.begin(); step != info::ecf->structural_mechanics_2d.load_steps_settings.end(); ++step) {
			_withFETI |= step->second.solver == LoadStepSolverConfiguration::SOLVER::FETI;
		}
		break;
	case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
		for (auto step = info::ecf->structural_mechanics_3d.load_steps_settings.begin(); step != info::ecf->structural_mechanics_3d.load_steps_settings.end(); ++step) {
			_withFETI |= step->second.solver == LoadStepSolverConfiguration::SOLVER::FETI;
		}
		break;
	default:
		eslog::globalerror("Mesh: Not implemented physics.\n");
		exit(0);
	}

	// check BEM
	for (auto it = info::ecf->getPhysics()->discretization.begin(); it != info::ecf->getPhysics()->discretization.end(); ++it) {
		if (it->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
			_withBEM = info::ecf->getPhysics()->dimension == DIMENSION::D3;
		}
	}

	// check whether we need to 'increase' boundary dimensions
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
				bregion(bc->first)->dimension = dimension;
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

		for (auto ls = ht->load_steps_settings.begin(); ls != ht->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.heat_flow.begin(); bc != ls->second.heat_flow.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
			for (auto bc = ls->second.heat_flux.begin(); bc != ls->second.heat_flux.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
			for (auto bc = ls->second.convection.begin(); bc != ls->second.convection.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
			for (auto bc = ls->second.diffuse_radiation.begin(); bc != ls->second.diffuse_radiation.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
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

		for (auto ls = te->load_steps_settings.begin(); ls != te->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.heat_transfer.heat_flow.begin(); bc != ls->second.heat_transfer.heat_flow.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
			for (auto bc = ls->second.heat_transfer.heat_flux.begin(); bc != ls->second.heat_transfer.heat_flux.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
			for (auto bc = ls->second.heat_transfer.convection.begin(); bc != ls->second.heat_transfer.convection.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
			for (auto bc = ls->second.heat_transfer.diffuse_radiation.begin(); bc != ls->second.heat_transfer.diffuse_radiation.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
			for (auto bc = ls->second.structural_mechanics.normal_pressure.begin(); bc != ls->second.structural_mechanics.normal_pressure.end(); ++bc) {
				bregion(bc->first)->dimension = dimension;
			}
		}
	}

	if (info::ecf->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D) {
		if (info::ecf->structural_mechanics_3d.discretization.size() && info::ecf->structural_mechanics_3d.discretization.begin()->second == PhysicsConfiguration::DISCRETIZATION::FEM_TDNNS) {
			_withEdgeDual = true;
		}
	}

	if (info::ecf->input_type == ECF::INPUT_TYPE::EXTERNAL_FILE) {
		_omitClusterization = false;
		_omitDecomposition = false;
	}

	if (info::ecf->input_type == ECF::INPUT_TYPE::GENERATOR) {
		_omitDecomposition = info::ecf->generator.uniform_clusters && info::ecf->generator.uniform_domains;
		_omitClusterization = info::ecf->generator.uniform_clusters;
	}

	if (_withBEM) {
		// TODO: BEM does not always need separate regions
		info::ecf->input.decomposition.separate_materials = true;
		info::ecf->input.decomposition.separate_regions = true;
	}
	if (
			info::ecf->input.decomposition.separate_materials ||
			info::ecf->input.decomposition.separate_regions ||
			info::ecf->input.decomposition.separate_etypes) {
		_omitDecomposition = false;
	}
}

void Mesh::preprocess()
{
	profiler::syncstart("meshing");

	analyze();

	eslog::startln("MESH: PREPROCESSING STARTED", "MESHING");
	setMaterials();
	eslog::checkpointln("MESH: MATERIALS FILLED");

	if (!_omitClusterization) {
		mesh::reclusterize();
		profiler::synccheckpoint("reclusterize");
	}

	esint msize;
	Communication::allReduce(&elements->size, &msize, 1, MPITools::getType<esint>().mpitype, MPI_MIN);
	if (msize == 0) {
		eslog::globalerror("ESPRESO quit computation: process without any elements detected.\n");
	}

	mesh::partitiate(preferedDomains, _omitDecomposition);

	profiler::synccheckpoint("domain_decomposition");
	eslog::checkpointln("MESH: MESH DECOMPOSED");


	for (size_t r = 0; r < boundaryRegions.size(); ++r) {
		if (boundaryRegions[r]->originalDimension < boundaryRegions[r]->dimension) {
			mesh::computeBoundaryElementsFromNodes(boundaryRegions[r]);
		}
	}
	eslog::checkpointln("MESH: BOUNDARY REGIONS COMPOSED");

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

	if (_withFETI) {
		mesh::computeNodeDomainDistribution();
		mesh::computeLocalIndices();
		profiler::synccheckpoint("preprocess_domains");
		eslog::checkpointln("MESH: ELEMENTS DOMAIN INDICES COMPUTED");
	}

	if (_withBEM) {
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

	if (_withFETI) {
		mesh::computeDomainDual();
		profiler::synccheckpoint("compute_domain_dual");
	}

	if (_withEdgeDual) {
		mesh::computeElementsEdgeNeighbors();
		profiler::synccheckpoint("compute_edge_neighbors");
		eslog::checkpointln("MESH: EDGE NEIGHBORS COMPUTED");
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
		packedSize += utils::packedSize(_omitClusterization);
		packedSize += utils::packedSize(_omitDecomposition);
		packedSize += utils::packedSize(_withGUI);
		packedSize += utils::packedSize(_withGUI);
		packedSize += utils::packedSize(_withFETI);
		packedSize += utils::packedSize(_withBEM);
		packedSize += utils::packedSize(_withEdgeDual);
	}

	Communication::broadcast(&packedSize, sizeof(size_t), MPI_BYTE, 0, MPITools::instances);
	char *buffer = new char[packedSize];

	if (info::mpi::irank == 0) {
		char *p = buffer;
		utils::pack(dimension, p);
		utils::pack(preferedDomains, p);

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
		utils::pack(_omitClusterization, p);
		utils::pack(_omitDecomposition, p);
		utils::pack(_withGUI, p);
		utils::pack(_withFETI, p);
		utils::pack(_withBEM, p);
		utils::pack(_withEdgeDual, p);
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
		utils::unpack(_omitClusterization, p);
		utils::unpack(_omitDecomposition, p);
		utils::unpack(_withGUI, p);
		utils::unpack(_withFETI, p);
		utils::unpack(_withBEM, p);
		utils::unpack(_withEdgeDual, p);

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
		switch (boundaryRegions[r]->originalDimension) {
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
			if (boundaryRegions[r]->originalDimension == 2) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), Parser::stringwithcommas(boundaryRegions[r]->totalsize).c_str(), Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "EDGES REGIONS SIZES", eregs, namesize - 26, " ", Parser::stringwithcommas(nedges).c_str(), Parser::stringwithcommas(nedgenodes).c_str());
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->originalDimension == 1) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), Parser::stringwithcommas(boundaryRegions[r]->totalsize).c_str(), Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "NODES REGIONS SIZES", nregs, namesize - 26, " ", " ", Parser::stringwithcommas(nnodes).c_str());
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->originalDimension == 0) {
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
			eslog::info("mesh: region=%s, dimension=%d, elements=%d, nodes=%d\n", boundaryRegions[r]->name.c_str(), boundaryRegions[r]->originalDimension, boundaryRegions[r]->totalsize, boundaryRegions[r]->nodeInfo.totalSize);
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

struct __meshinfo__ {
	struct op_min { template<typename T> T operator()(const T&i, const T&j) const { return std::min(i, j); } };
	struct op_max { template<typename T> T operator()(const T&i, const T&j) const { return std::max(i, j); } };
	struct op_sum { template<typename T> T operator()(const T&i, const T&j) const { return i + j; } };

	template <typename T>
	struct domain_t {
		T elements = 0, uniquenodes = 0, nodes = 0, neighbors = 0, ninterface = 0, einterface = 0;
		double nratio = 0, eratio = 0;

		void set(const T &value)
		{
			elements = uniquenodes = nodes = neighbors = ninterface = einterface = value;
			nratio = eratio = value;
		}

		template<typename TT>
		domain_t& operator=(const domain_t<TT> &other)
		{
			elements = other.elements;
			uniquenodes = other.uniquenodes;
			nodes = other.nodes;
			neighbors = other.neighbors;
			ninterface = other.ninterface;
			einterface = other.einterface;
			nratio = other.nratio;
			eratio = other.eratio;
			return *this;
		}

		template <typename OP>
		void apply(const domain_t &other, const OP &op)
		{
			elements = op(elements, other.elements);
			uniquenodes = op(uniquenodes, other.uniquenodes);
			nodes = op(nodes, other.nodes);
			neighbors = op(neighbors, other.neighbors);
			ninterface = op(ninterface, other.ninterface);
			einterface = op(einterface, other.einterface);
			nratio = op(nratio, other.nratio);
			eratio = op(eratio, other.eratio);
		}

		void avg(const domain_t<esint> &sum)
		{
			elements /= sum.elements;
			uniquenodes /= sum.uniquenodes;
			nodes /= sum.nodes;
			neighbors /= sum.neighbors;
			ninterface /= sum.ninterface;
			einterface /= sum.einterface;
			nratio /= sum.nratio;
			eratio /= sum.eratio;
		}

		void var(const domain_t<esint> &sum)
		{
			if (sum.elements - 1 > 0) elements = std::sqrt(elements * (1. / (sum.elements - 1)));
			if (sum.uniquenodes - 1 > 0) uniquenodes = std::sqrt(uniquenodes * (1. / (sum.uniquenodes - 1)));
			if (sum.nodes - 1 > 0) nodes = std::sqrt(nodes * (1. / (sum.nodes - 1)));
			if (sum.neighbors - 1 > 0) neighbors = std::sqrt(neighbors * (1. / (sum.neighbors - 1)));
			if (sum.ninterface - 1 > 0) ninterface = std::sqrt(ninterface * (1. / (sum.ninterface - 1)));
			if (sum.einterface - 1 > 0) einterface = std::sqrt(einterface * (1. / (sum.einterface - 1)));
			if (sum.nratio - 1 > 0) nratio = std::sqrt(nratio * (1. / (sum.nratio - 1)));
			if (sum.eratio - 1 > 0) eratio = std::sqrt(eratio * (1. / (sum.eratio - 1)));
		}
	};
	template <typename T>
	struct cluster_t: public domain_t<T> {
		T domains = 0;

		void set(const T &value)
		{
			domain_t<T>::set(value);
			domains = value;
		}

		template<typename TT>
		cluster_t& operator=(const cluster_t<TT> &other)
		{
			domain_t<T>::operator=(other);
			domains = other.domains;
			return *this;
		}

		template<typename TT>
		cluster_t& operator=(const domain_t<TT> &other)
		{
			domain_t<T>::operator=(other);
			return *this;
		}

		template <typename OP>
		void apply(const cluster_t &other, const OP &op)
		{
			domain_t<T>::apply(other, op);
			domains = op(domains, other.domains);
		}

		void avg(const cluster_t<esint> &sum)
		{
			domain_t<T>::avg(sum);
			domains /= sum.domains;
		}

		void var(const cluster_t<esint> &sum)
		{
			domain_t<T>::var(sum);
			if (sum.domains - 1 > 0) domains = std::sqrt(domains * (1. / (sum.domains - 1)));
		}
	};
	template <typename T>
	struct mpi_t: public cluster_t<T> {
		T clusters = 0;

		void set(const T &value)
		{
			cluster_t<T>::set(value);
			clusters = value;
		}

		template<typename TT>
		mpi_t& operator=(const mpi_t<TT> &other)
		{
			cluster_t<T>::operator=(other);
			clusters = other.clusters;
			return *this;
		}

		template <typename OP>
		void apply(const mpi_t &other, const OP &op)
		{
			cluster_t<T>::apply(other, op);
			clusters = op(clusters, other.clusters);
		}

		void avg(const mpi_t<esint> &sum)
		{
			cluster_t<T>::avg(sum);
			clusters /= sum.clusters;
		}

		void var(const mpi_t<esint> &sum)
		{
			cluster_t<T>::var(sum);
			if (sum.clusters - 1 > 0) clusters = std::sqrt(clusters * (1. / (sum.clusters - 1)));
		}
	};

	template<template<typename> typename C, typename T>
	struct stats_t {
		C<T> min, max, sum;
		C<double> avg;

		template<template<typename> typename CC, typename TT>
		stats_t<C, T>& operator=(const stats_t<CC, TT> &other)
		{
			min = other.min;
			max = other.max;
			sum = other.sum;
			avg = other.avg;
			return *this;
		}

		void set(const C<T> &value)
		{
			avg = min = max = sum = value;
		}

		void apply(const stats_t<C, T> &other)
		{
			min.apply(other.min, op_min());
			max.apply(other.max, op_max());
			avg.apply(other.avg, op_sum());
			sum.apply(other.sum, op_sum());
		}
	};

	struct value_t {
		mpi_t<esint> mpi;
		cluster_t<esint> cluster;
		domain_t<esint> domain;
	} value;

	struct variance_t {
		mpi_t<double> mpi;
		cluster_t<double> cluster;
		domain_t<double> domain;
	} variance;

	struct base_t {
		stats_t<mpi_t, esint> mpi;
		stats_t<cluster_t, esint> cluster;
		stats_t<domain_t, esint> domain;

		void apply(const base_t &other)
		{
			mpi.apply(other.mpi);
			cluster.apply(other.cluster);
			domain.apply(other.domain);
		}
	} stats;

	void avg()
	{
		stats.mpi.avg.avg(stats.mpi.sum);
		stats.cluster.avg.avg(stats.cluster.sum);
		stats.domain.avg.avg(stats.domain.sum);
	}

	void var()
	{
		variance.mpi.var(stats.mpi.sum);
		variance.cluster.var(stats.cluster.sum);
		variance.domain.var(stats.domain.sum);
	}

	MPI_Datatype base_mpi_t;
	MPI_Op base_mpi_op;

	static void op(void *in, void *out, int *len, MPI_Datatype *datatype)
	{
		for (int i = 0; i < *len; i++) {
			(static_cast<base_t*>(out) + i)->apply(*(static_cast<base_t*>(in) + i));
		}
	}

	__meshinfo__()
	{
		MPI_Type_contiguous(sizeof(stats), MPI_BYTE, &base_mpi_t);
		MPI_Type_commit(&base_mpi_t);
		MPI_Op_create(op, true, &base_mpi_op);
	}

	~__meshinfo__()
	{
		MPI_Type_free(&base_mpi_t);
		MPI_Op_free(&base_mpi_op);
	}
};

void Mesh::printDecompositionStatistics()
{
	__meshinfo__ mesh;
	std::vector<esint> dnodes(elements->ndomains), dninterface(elements->ndomains), delements(elements->ndomains), deinterface(elements->ndomains), dneighs(elements->ndomains);
	std::vector<double> dnratio(elements->ndomains), deratio(elements->ndomains);
	std::vector<esint> cnodes(elements->nclusters), cninterface(elements->nclusters), celements(elements->nclusters), ceinterface(elements->nclusters), cneighs(elements->nclusters), cdomains(elements->nclusters);
	std::vector<double> cnratio(elements->nclusters), ceratio(elements->nclusters);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		esint inodes = 0;
		for (auto ranks = nodes->ranks->begin(t); ranks != nodes->ranks->end(t); ++ranks) {
			if (ranks->size() > 1) {
				++inodes;
			}
		}
		#pragma omp atomic
		mesh.value.mpi.ninterface += inodes;

		esint ielements = 0;
		for (auto neighs = elements->faceNeighbors->begin(t); neighs != elements->faceNeighbors->end(t); ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1 && (*n < elements->offset || elements->offset + elements->size <= *n)) {
					++ielements;
				}
			}
		}
		#pragma omp atomic
		mesh.value.mpi.einterface += ielements;

		if (_withFETI) {
			std::vector<esint> tdnodes(elements->ndomains), tdninterface(elements->ndomains);
			for (auto dmap = nodes->domains->begin(t); dmap != nodes->domains->end(t); ++dmap) {
				if (dmap->size() == 1) {
					++tdnodes[dmap->front() - elements->firstDomain];
				} else {
					for (auto d = dmap->begin(); d != dmap->end(); ++d) {
						if (elements->firstDomain <= *d && *d < elements->firstDomain + elements->ndomains) {
							++tdnodes[*d - elements->firstDomain];
							++tdninterface[*d - elements->firstDomain];
						}
					}
				}
			}

			std::vector<esint> tdeinterface(elements->ndomains);
			for (esint d = elements->domainDistribution[t]; d < elements->domainDistribution[t + 1]; ++d) {
				auto begin = elements->faceNeighbors->begin() + elements->elementsDistribution[d];
				auto end = elements->faceNeighbors->begin() + elements->elementsDistribution[d + 1];
				for (auto neighs = begin; neighs != end; ++neighs) {
					esint size = 0;
					for (auto n = neighs->begin(); n != neighs->end(); ++n) {
						if (*n != -1 && (*n < elements->offset + elements->elementsDistribution[d] || elements->offset + elements->elementsDistribution[d + 1] <= *n)) {
							++tdeinterface[d];
							++size;
						}
					}
				}
			}

			for (esint d = 0; d < elements->ndomains; ++d) {
				dnodes[d] += tdnodes[d];
				dninterface[d] += tdninterface[d];
				deinterface[d] += tdeinterface[d];
			}
		}
	}

	mesh.value.mpi.clusters    = elements->nclusters;
	mesh.value.mpi.domains     = elements->ndomains;
	mesh.value.mpi.elements    = elements->size;
	mesh.value.mpi.uniquenodes = nodes->uniqInfo.size;
	mesh.value.mpi.nodes       = nodes->size;
	mesh.value.mpi.nratio      = (double)mesh.value.mpi.ninterface / mesh.value.mpi.nodes;
	mesh.value.mpi.eratio      = (double)mesh.value.mpi.einterface / mesh.value.mpi.elements;
	mesh.value.mpi.neighbors   = neighbors.size();

	mesh.stats.mpi.set(mesh.value.mpi);
	mesh.stats.mpi.sum.set(1);

	if (_withFETI) {
		mesh.value.domain.elements += elements->size;
		for (esint d = 0; d < elements->ndomains; ++d) {
			mesh.value.domain.nodes += dnodes[d];
			mesh.value.domain.ninterface += dninterface[d];
			dnratio[d] = (double)dninterface[d] / dnodes[d];
			mesh.value.domain.nratio += dnratio[d];
			delements[d] = elements->elementsDistribution[d + 1] - elements->elementsDistribution[d];
			mesh.value.domain.einterface += deinterface[d];
			deratio[d] = (double)deinterface[d] / delements[d];
			mesh.value.domain.eratio += deratio[d];
			dneighs[d] = (FETIData->fullDomainDual->begin() + d)->size();
			mesh.value.domain.neighbors += dneighs[d];
		}
		mesh.stats.domain.set(mesh.value.domain);
		{
			auto minmax = std::minmax_element(dnodes.begin(), dnodes.end());
			mesh.stats.domain.min.nodes = *minmax.first;
			mesh.stats.domain.max.nodes = *minmax.second;
		}
		{
			auto minmax = std::minmax_element(dninterface.begin(), dninterface.end());
			mesh.stats.domain.min.ninterface = *minmax.first;
			mesh.stats.domain.max.ninterface = *minmax.second;
		}
		{
			auto minmax = std::minmax_element(dnratio.begin(), dnratio.end());
			mesh.stats.domain.min.nratio = *minmax.first;
			mesh.stats.domain.max.nratio = *minmax.second;
		}
		{
			auto minmax = std::minmax_element(delements.begin(), delements.end());
			mesh.stats.domain.min.elements = *minmax.first;
			mesh.stats.domain.max.elements = *minmax.second;
		}
		{
			auto minmax = std::minmax_element(deinterface.begin(), deinterface.end());
			mesh.stats.domain.min.einterface = *minmax.first;
			mesh.stats.domain.max.einterface = *minmax.second;
		}
		{
			auto minmax = std::minmax_element(deratio.begin(), deratio.end());
			mesh.stats.domain.min.eratio = *minmax.first;
			mesh.stats.domain.max.eratio = *minmax.second;
		}
		{
			auto minmax = std::minmax_element(dneighs.begin(), dneighs.end());
			mesh.stats.domain.min.neighbors = *minmax.first;
			mesh.stats.domain.max.neighbors = *minmax.second;
		}
		mesh.stats.domain.sum.set(elements->ndomains);
		mesh.value.cluster = mesh.value.mpi;
		mesh.stats.cluster = mesh.stats.mpi;
		if (elements->nclusters > 1) {
			mesh.stats.cluster.min = mesh.stats.domain.max;
			mesh.stats.cluster.max = mesh.stats.domain.min;
			for (esint d = 0; d < elements->ndomains; ++d) {
				++cdomains[elements->clusters[d]];
				celements[elements->clusters[d]] += delements[d];
				cnodes[elements->clusters[d]] += dnodes[d];
				cninterface[elements->clusters[d]] += dninterface[d];
				cnratio[elements->clusters[d]] += dnratio[d];
				ceinterface[elements->clusters[d]] += deinterface[d];
				ceratio[elements->clusters[d]] += deratio[d];
				cneighs[elements->clusters[d]] += dneighs[d];
			}
			{
				auto minmax = std::minmax_element(cdomains.begin(), cdomains.end());
				mesh.stats.cluster.min.domains = *minmax.first;
				mesh.stats.cluster.max.domains = *minmax.second;
			}
			{
				auto minmax = std::minmax_element(cnodes.begin(), cnodes.end());
				mesh.stats.cluster.min.nodes = *minmax.first;
				mesh.stats.cluster.max.nodes = *minmax.second;
			}
			{
				auto minmax = std::minmax_element(cninterface.begin(), cninterface.end());
				mesh.stats.cluster.min.ninterface = *minmax.first;
				mesh.stats.cluster.max.ninterface = *minmax.second;
			}
			{
				auto minmax = std::minmax_element(cnratio.begin(), cnratio.end());
				mesh.stats.cluster.min.nratio = *minmax.first;
				mesh.stats.cluster.max.nratio = *minmax.second;
			}
			{
				auto minmax = std::minmax_element(celements.begin(), celements.end());
				mesh.stats.cluster.min.elements = *minmax.first;
				mesh.stats.cluster.max.elements = *minmax.second;
			}
			{
				auto minmax = std::minmax_element(ceinterface.begin(), ceinterface.end());
				mesh.stats.cluster.min.einterface = *minmax.first;
				mesh.stats.cluster.max.einterface = *minmax.second;
			}
			{
				auto minmax = std::minmax_element(ceratio.begin(), ceratio.end());
				mesh.stats.cluster.min.eratio = *minmax.first;
				mesh.stats.cluster.max.eratio = *minmax.second;
			}
			{
				auto minmax = std::minmax_element(cneighs.begin(), cneighs.end());
				mesh.stats.cluster.min.neighbors = *minmax.first;
				mesh.stats.cluster.max.neighbors = *minmax.second;
			}
			mesh.stats.cluster.sum.set(elements->nclusters);
		}
	}

	Communication::allReduce(&mesh.stats, NULL, 1, mesh.base_mpi_t, mesh.base_mpi_op);
	mesh.avg();
	mesh.variance.mpi.elements    = std::pow(mesh.stats.mpi.avg.elements    - mesh.value.mpi.elements, 2);
	mesh.variance.mpi.uniquenodes = std::pow(mesh.stats.mpi.avg.uniquenodes - mesh.value.mpi.uniquenodes, 2);
	mesh.variance.mpi.nodes       = std::pow(mesh.stats.mpi.avg.nodes       - mesh.value.mpi.nodes, 2);
	mesh.variance.mpi.neighbors   = std::pow(mesh.stats.mpi.avg.neighbors   - mesh.value.mpi.neighbors, 2);
	mesh.variance.mpi.ninterface  = std::pow(mesh.stats.mpi.avg.ninterface  - mesh.value.mpi.ninterface, 2);
	mesh.variance.mpi.einterface  = std::pow(mesh.stats.mpi.avg.einterface  - mesh.value.mpi.einterface, 2);
	mesh.variance.mpi.nratio      = std::pow(mesh.stats.mpi.avg.nratio      - mesh.value.mpi.nratio, 2);
	mesh.variance.mpi.eratio      = std::pow(mesh.stats.mpi.avg.eratio      - mesh.value.mpi.eratio, 2);
	mesh.variance.mpi.domains     = std::pow(mesh.stats.mpi.avg.domains     - mesh.value.mpi.domains, 2);
	mesh.variance.mpi.clusters    = std::pow(mesh.stats.mpi.avg.clusters    - mesh.value.mpi.clusters, 2);
	if (mesh.stats.mpi.max.clusters == info::mpi::size) {
		mesh.variance.cluster = mesh.variance.mpi;
	} else {
		for (esint c = 0; c < elements->nclusters; ++c) {
			mesh.variance.cluster.elements    = std::pow(mesh.stats.cluster.avg.elements    - celements[c], 2);
			mesh.variance.cluster.nodes       = std::pow(mesh.stats.cluster.avg.nodes       - cnodes[c], 2);
			mesh.variance.cluster.ninterface  = std::pow(mesh.stats.cluster.avg.ninterface  - cninterface[c], 2);
			mesh.variance.cluster.einterface  = std::pow(mesh.stats.cluster.avg.einterface  - ceinterface[c], 2);
			mesh.variance.cluster.nratio      = std::pow(mesh.stats.cluster.avg.nratio      - cnratio[c], 2);
			mesh.variance.cluster.eratio      = std::pow(mesh.stats.cluster.avg.eratio      - ceratio[c], 2);
			mesh.variance.cluster.neighbors   = std::pow(mesh.stats.cluster.avg.neighbors   - cneighs[c], 2);
		}
	}
	for (esint d = 0; d < elements->ndomains; ++d) {
		mesh.variance.domain.elements    = std::pow(mesh.stats.domain.avg.elements    - delements[d], 2);
		mesh.variance.domain.nodes       = std::pow(mesh.stats.domain.avg.nodes       - dnodes[d], 2);
		mesh.variance.domain.ninterface  = std::pow(mesh.stats.domain.avg.ninterface  - dninterface[d], 2);
		mesh.variance.domain.einterface  = std::pow(mesh.stats.domain.avg.einterface  - deinterface[d], 2);
		mesh.variance.domain.nratio      = std::pow(mesh.stats.domain.avg.nratio      - dnratio[d], 2);
		mesh.variance.domain.eratio      = std::pow(mesh.stats.domain.avg.eratio      - deratio[d], 2);
		mesh.variance.domain.neighbors   = std::pow(mesh.stats.domain.avg.neighbors   - dneighs[d], 2);
	}

	Communication::allReduce(&mesh.variance, NULL, sizeof(mesh.variance) / sizeof(double), MPI_DOUBLE, MPI_SUM);
	mesh.var();

	eslog::info(" ================================== DECOMPOSITION STATISTICS ================================= \n");
	if (!_withFETI || mesh.stats.mpi.max.clusters > 1) {
		eslog::info(" ============================================================================================= \n");
		eslog::info("  PER MPI STATISTICS         %12s %12s %12s %12s %12s \n", "MIN", "MAX", "AVG", "VAR", "IMBALANCE");
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		eslog::info("  CLUSTERS                 : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.clusters,
				mesh.stats.mpi.max.clusters,
				mesh.stats.mpi.avg.clusters,
				mesh.variance.mpi.clusters,
				mesh.stats.mpi.min.clusters ? (double)mesh.stats.mpi.max.clusters / mesh.stats.mpi.min.clusters : 0.0);
		eslog::info("  DOMAINS                  : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.domains,
				mesh.stats.mpi.max.domains,
				mesh.stats.mpi.avg.domains,
				mesh.variance.mpi.domains,
				mesh.stats.mpi.min.domains ? (double)mesh.stats.mpi.max.domains / mesh.stats.mpi.min.domains : 0.0);
		eslog::info("  UNIQUE NODES             : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.uniquenodes).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.uniquenodes).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.uniquenodes).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.uniquenodes).c_str(),
				mesh.stats.mpi.min.uniquenodes ? (double)mesh.stats.mpi.max.uniquenodes / mesh.stats.mpi.min.uniquenodes : 0.0);
		eslog::info("  NODES                    : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.nodes).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.nodes).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.nodes).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.nodes).c_str(),
				mesh.stats.mpi.min.nodes ? (double)mesh.stats.mpi.max.nodes / mesh.stats.mpi.min.nodes : 0.0);
		eslog::info("  INTERFACE NODES          : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.ninterface).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.ninterface).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.ninterface).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.ninterface).c_str(),
				mesh.stats.mpi.min.ninterface ? (double)mesh.stats.mpi.max.ninterface / mesh.stats.mpi.min.ninterface : 0.0);
		eslog::info("  INTERFACE NODES RATIO    : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.nratio,
				mesh.stats.mpi.max.nratio,
				mesh.stats.mpi.avg.nratio,
				mesh.variance.mpi.nratio,
				mesh.stats.mpi.min.nratio ? mesh.stats.mpi.max.nratio / mesh.stats.mpi.min.nratio : 0.0);
		eslog::info("  ELEMENTS                 : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.elements).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.elements).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.elements).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.elements).c_str(),
				mesh.stats.mpi.min.elements ? (double)mesh.stats.mpi.max.elements / mesh.stats.mpi.min.elements : 0.0);
		eslog::info("  INTERFACE ELEMENTS       : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.einterface).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.einterface).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.einterface).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.einterface).c_str(),
				mesh.stats.mpi.min.einterface ? (double)mesh.stats.mpi.max.einterface / mesh.stats.mpi.min.einterface : 0);
		eslog::info("  INTERFACE ELEMENTS RATIO : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.eratio,
				mesh.stats.mpi.max.eratio,
				mesh.stats.mpi.avg.eratio,
				mesh.variance.mpi.eratio,
				mesh.stats.mpi.min.eratio ? mesh.stats.mpi.max.eratio / mesh.stats.mpi.min.eratio : 0.0);
		eslog::info("  NEIGHBORS                : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.neighbors,
				mesh.stats.mpi.max.neighbors,
				mesh.stats.mpi.avg.neighbors,
				mesh.variance.mpi.neighbors,
				mesh.stats.mpi.min.neighbors ? (double)mesh.stats.mpi.max.neighbors / mesh.stats.mpi.min.neighbors : 0.0);
	}

	if (_withFETI) {
		eslog::info(" ============================================================================================= \n");
		eslog::info("  PER CLUSTER STATISTICS     %12s %12s %12s %12s %12s \n", "MIN", "MAX", "AVG", "VAR", "IMBALANCE");
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		eslog::info("  DOMAINS                  : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.cluster.min.domains,
				mesh.stats.cluster.max.domains,
				mesh.stats.cluster.avg.domains,
				mesh.variance.cluster.domains,
				mesh.stats.cluster.min.domains ? (double)mesh.stats.cluster.max.domains / mesh.stats.cluster.min.domains : 0.0);
		if (mesh.stats.mpi.max.clusters == 1) {
			eslog::info("  UNIQUE NODES             : %12s %12s %12s %12s %12.2f \n",
					Parser::stringwithcommas(mesh.stats.cluster.min.uniquenodes).c_str(),
					Parser::stringwithcommas(mesh.stats.cluster.max.uniquenodes).c_str(),
					Parser::stringwithcommas((int)mesh.stats.cluster.avg.uniquenodes).c_str(),
					Parser::stringwithcommas((int)mesh.variance.cluster.uniquenodes).c_str(),
					mesh.stats.cluster.min.uniquenodes ? (double)mesh.stats.cluster.max.uniquenodes / mesh.stats.cluster.min.uniquenodes : 0.0);
			eslog::info("  NODES                    : %12s %12s %12s %12s %12.2f \n",
					Parser::stringwithcommas(mesh.stats.cluster.min.nodes).c_str(),
					Parser::stringwithcommas(mesh.stats.cluster.max.nodes).c_str(),
					Parser::stringwithcommas((int)mesh.stats.cluster.avg.nodes).c_str(),
					Parser::stringwithcommas((int)mesh.variance.cluster.nodes).c_str(),
					mesh.stats.cluster.min.nodes ? (double)mesh.stats.cluster.max.nodes / mesh.stats.cluster.min.nodes : 0.0);
			eslog::info("  INTERFACE NODES          : %12s %12s %12s %12s %12.2f \n",
					Parser::stringwithcommas(mesh.stats.cluster.min.ninterface).c_str(),
					Parser::stringwithcommas(mesh.stats.cluster.max.ninterface).c_str(),
					Parser::stringwithcommas((int)mesh.stats.cluster.avg.ninterface).c_str(),
					Parser::stringwithcommas((int)mesh.variance.cluster.ninterface).c_str(),
					mesh.stats.cluster.min.ninterface ? (double)mesh.stats.cluster.max.ninterface / mesh.stats.cluster.min.ninterface : 0.0);
			eslog::info("  INTERFACE NODES RATIO    : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
					mesh.stats.cluster.min.nratio,
					mesh.stats.cluster.max.nratio,
					mesh.stats.cluster.avg.nratio,
					mesh.variance.cluster.nratio,
					mesh.stats.cluster.min.nratio ? mesh.stats.cluster.max.nratio / mesh.stats.cluster.min.nratio : 0.0);
		}
		eslog::info("  ELEMENTS                 : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.cluster.min.elements).c_str(),
				Parser::stringwithcommas(mesh.stats.cluster.max.elements).c_str(),
				Parser::stringwithcommas((int)mesh.stats.cluster.avg.elements).c_str(),
				Parser::stringwithcommas((int)mesh.variance.cluster.elements).c_str(),
				mesh.stats.cluster.min.elements ? (double)mesh.stats.cluster.max.elements / mesh.stats.cluster.min.elements : 0.0);
		if (mesh.stats.mpi.max.clusters == 1) {
			eslog::info("  INTERFACE ELEMENTS       : %12s %12s %12s %12s %12.2f \n",
					Parser::stringwithcommas(mesh.stats.cluster.min.einterface).c_str(),
					Parser::stringwithcommas(mesh.stats.cluster.max.einterface).c_str(),
					Parser::stringwithcommas((int)mesh.stats.cluster.avg.einterface).c_str(),
					Parser::stringwithcommas((int)mesh.variance.cluster.einterface).c_str(),
					mesh.stats.cluster.min.einterface ? (double)mesh.stats.cluster.max.einterface / mesh.stats.cluster.min.einterface : 0.0);
			eslog::info("  INTERFACE ELEMENTS RATIO : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
					mesh.stats.cluster.min.eratio,
					mesh.stats.cluster.max.eratio,
					mesh.stats.cluster.avg.eratio,
					mesh.variance.cluster.eratio,
					mesh.stats.cluster.min.eratio ? mesh.stats.cluster.max.eratio / mesh.stats.cluster.min.eratio : 0.0);
			if (mesh.stats.mpi.max.clusters == 1) {
				eslog::info("  NEIGHBORS                : %12d %12d %12.2f %12.2f %12.2f \n",
						mesh.stats.cluster.min.neighbors,
						mesh.stats.cluster.max.neighbors,
						mesh.stats.cluster.avg.neighbors,
						mesh.variance.cluster.neighbors,
						mesh.stats.cluster.min.neighbors ? (double)mesh.stats.cluster.max.neighbors / mesh.stats.cluster.min.neighbors : 0.0);
			}
		}

		eslog::info(" ============================================================================================= \n");
		eslog::info("  PER DOMAIN STATISTICS      %12s %12s %12s %12s %12s \n", "MIN", "MAX", "AVG", "VAR", "IMBALANCE");
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		eslog::info("  NODES                    : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.domain.min.nodes).c_str(),
				Parser::stringwithcommas(mesh.stats.domain.max.nodes).c_str(),
				Parser::stringwithcommas((int)mesh.stats.domain.avg.nodes).c_str(),
				Parser::stringwithcommas((int)mesh.variance.domain.nodes).c_str(),
				mesh.stats.domain.min.nodes ? (double)mesh.stats.domain.max.nodes / mesh.stats.domain.min.nodes : 0.0);
		eslog::info("  INTERFACE NODES          : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.domain.min.ninterface).c_str(),
				Parser::stringwithcommas(mesh.stats.domain.max.ninterface).c_str(),
				Parser::stringwithcommas((int)mesh.stats.domain.avg.ninterface).c_str(),
				Parser::stringwithcommas((int)mesh.variance.domain.ninterface).c_str(),
				mesh.stats.domain.min.ninterface ? (double)mesh.stats.domain.max.ninterface / mesh.stats.domain.min.ninterface : 0.0);
		eslog::info("  INTERFACE NODES RATIO    : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
				mesh.stats.domain.min.nratio,
				mesh.stats.domain.max.nratio,
				mesh.stats.domain.avg.nratio,
				mesh.variance.domain.nratio,
				mesh.stats.domain.min.nratio ? mesh.stats.domain.max.nratio / mesh.stats.domain.min.nratio : 0.0);
		eslog::info("  ELEMENTS                 : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.domain.min.elements).c_str(),
				Parser::stringwithcommas(mesh.stats.domain.max.elements).c_str(),
				Parser::stringwithcommas((int)mesh.stats.domain.avg.elements).c_str(),
				Parser::stringwithcommas((int)mesh.variance.domain.elements).c_str(),
				mesh.stats.domain.min.elements ? (double)mesh.stats.domain.max.elements / mesh.stats.domain.min.elements : 0.0);
		eslog::info("  INTERFACE ELEMENTS       : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.domain.min.einterface).c_str(),
				Parser::stringwithcommas(mesh.stats.domain.max.einterface).c_str(),
				Parser::stringwithcommas((int)mesh.stats.domain.avg.einterface).c_str(),
				Parser::stringwithcommas((int)mesh.variance.domain.einterface).c_str(),
				mesh.stats.domain.min.einterface ? (double)mesh.stats.domain.max.einterface / mesh.stats.domain.min.einterface : 0.0);
		eslog::info("  INTERFACE ELEMENTS RATIO : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
				mesh.stats.domain.min.eratio,
				mesh.stats.domain.max.eratio,
				mesh.stats.domain.avg.eratio,
				mesh.variance.domain.eratio,
				mesh.stats.domain.min.eratio ? mesh.stats.domain.max.eratio / mesh.stats.domain.min.eratio : 0.0);
		eslog::info("  NEIGHBORS DOMAINS        : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.domain.min.neighbors,
				mesh.stats.domain.max.neighbors,
				mesh.stats.domain.avg.neighbors,
				mesh.variance.domain.neighbors,
				mesh.stats.domain.min.neighbors ? (double)mesh.stats.domain.max.neighbors / mesh.stats.domain.min.neighbors : 0.0);
	}
	eslog::info(" ============================================================================================= \n");

	if (info::ecf->output.store_decomposition > 1) {
		NamedData *data;
		if (!_withFETI || mesh.stats.mpi.max.clusters > 1) {
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_CLUSTERS");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.clusters);
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_NODES");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.nodes);
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_INODES");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.ninterface);
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_INODES_RATIO");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.nratio);
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_ELEMENTS");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.elements);
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_IELEMENTS");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.einterface);
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_IELEMENTS_RATIO");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.eratio);
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_NEIGHBORS");
			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.neighbors);
		}
		if (_withFETI) {
			if (mesh.stats.mpi.max.clusters == 1) {
				data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "CLUSTER_NODES");
				std::fill(data->data.begin(), data->data.end(), mesh.value.cluster.nodes);
				data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "CLUSTER_INODES");
				std::fill(data->data.begin(), data->data.end(), mesh.value.cluster.ninterface);
				data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "CLUSTER_INODES_RATIO");
				std::fill(data->data.begin(), data->data.end(), mesh.value.cluster.nratio);
			}
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "CLUSTER_ELEMENTS");
			if (elements->nclusters == 1) {
				std::fill(data->data.begin(), data->data.end(), mesh.value.cluster.elements);
			} else {
				for (esint d = 0; d < elements->ndomains; ++d) {
					std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], celements[elements->clusters[d]]);
				}
			}
			if (mesh.stats.mpi.max.clusters == 1) {
				data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "CLUSTER_IELEMENTS");
				std::fill(data->data.begin(), data->data.end(), mesh.value.cluster.einterface);
				data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "CLUSTER_IELEMENTS_RATIO");
				std::fill(data->data.begin(), data->data.end(), mesh.value.cluster.eratio);
				data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "CLUSTER_NEIGHBORS");
				std::fill(data->data.begin(), data->data.end(), mesh.value.cluster.neighbors);
			}

			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_NODES");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], dnodes[d]);
			}
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_INODES");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], dninterface[d]);
			}
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_INODES_RATIO");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], dnratio[d]);
			}
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_ELEMENTS");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], delements[d]);
			}
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_IELEMENTS");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], deinterface[d]);
			}
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_IELEMENTS_RATIO");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], deratio[d]);
			}

			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_IELEMENTS_RATIO");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], deratio[d]);
			}
			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "DOMAIN_NEIGHBORS");
			for (esint d = 0; d < elements->ndomains; ++d) {
				std::fill(data->data.begin() + elements->elementsDistribution[d], data->data.begin() + elements->elementsDistribution[d + 1], dneighs[d]);
			}
		}
	}
}

