
#include "heattransfer.kernel.hpp"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/elementsregionstore.h"

#include "basis/utilities/print.h"

using namespace espreso;

void HeatTransferKernel::nextSubstep()
{
	computeCoodinates();
	initTemperature();
	if (info::mesh->dimension == 2) {
		computeThickness();
	}
	computeJacobian();
	computeConductivity();
//	computeAdvection();

	computeStiffness();
}

void HeatTransferKernel::solutionChanged()
{

}

bool HeatTransferKernel::boundaryWithSettings(size_t rindex)
{
	return false;
}

template<int code, int nodes, int gpsize>
void HeatTransferKernel::isotropicK(const Builder &builder, InstanceFiller &filler)
{
//	filler.DOFs = nodes;
//	filler.Ke.resize(nodes, nodes);
//
//	filler.insertK = true;
//	filler.insertM = false;
//	filler.insertC = false;
//	filler.insertR = false;
//	filler.insertF = false;
//
//	const std::vector<double> &weighFactor = *Mesh::edata[code].weighFactor;
//
//	const double * K = (_isotropicK.data->begin() + filler.interval)->data();
//	int Kinc = _isotropicK.isconstant[filler.interval] ? 0 : 1;
//	const double *det = _detJ.data->datatarray().data();
//	const double *dND = _dND.data->datatarray().data();
//
//	if (info::mesh->dimension == 2) {
//		const double * thickness = (_thickness.data->begin() + filler.interval)->data();
//		int thicknessinc = _thickness.isconstant[filler.interval] ? 0 : 1;
//		for (esint e = filler.begin; e < filler.end; ++e) {
//			filler.Ke.fill(0);
//
//			for (int gp = 0; gp < gpsize; ++gp, ++det, dND += 2 * nodes, K += Kinc, thickness += thicknessinc) {
//				KMN2M2N<nodes>(*K * *thickness * *det * weighFactor[gp], dND, filler.Ke.vals);
//			}
//
//			filler.insert();
//		}
//	}
//	if (info::mesh->dimension == 3) {
//		for (esint e = filler.begin; e < filler.end; ++e) {
//			filler.Ke.fill(0);
//
//			for (int gp = 0; gp < gpsize; ++gp, ++det, dND += 3 * nodes, K += Kinc) {
//				KMN3M3N<nodes>(*K * *det * weighFactor[gp], dND, filler.Ke.vals);
//			}
//
//			filler.insert();
//		}
//	}
}

template<int code, int nodes, int gpsize>
void HeatTransferKernel::generalK(const Builder &builder, InstanceFiller &filler)
{
//	filler.DOFs = nodes;
//	filler.Ke.resize(nodes, nodes);
//
//	filler.insertK = true;
//	filler.insertM = false;
//	filler.insertC = false;
//	filler.insertR = false;
//	filler.insertF = false;
//
//	const std::vector<double> &weighFactor = *Mesh::edata[code].weighFactor;
//
//	const double * K = (_K.data->begin() + filler.interval)->data();
//	int Kinc = _K.isconstant[filler.interval] ? 0 : info::mesh->dimension * info::mesh->dimension;
//	const double *det = _detJ.data->datatarray().data();
//	const double *dND = _dND.data->datatarray().data();
//
//	if (info::mesh->dimension == 2) {
//		const double * thickness = (_thickness.data->begin() + filler.interval)->data();
//		int thicknessinc = _thickness.isconstant[filler.interval] ? 0 : 1;
//		for (esint e = filler.begin; e < filler.end; ++e) {
//			filler.Ke.fill(0);
//
//			for (int gp = 0; gp < gpsize; ++gp, ++det, dND += 2 * nodes, K += Kinc, thickness += thicknessinc) {
//				KMN2M22M2N<nodes>(*thickness * *det * weighFactor[gp], K, dND, filler.Ke.vals);
//			}
//
//			filler.insert();
//		}
//	}
//
//	if (info::mesh->dimension == 3) {
//		for (esint e = filler.begin; e < filler.end; ++e) {
//			filler.Ke.fill(0);
//
//			for (int gp = 0; gp < gpsize; ++gp, ++det, dND += 3 * nodes, K += Kinc) {
//				KMN3M33M3N<nodes>(*det * weighFactor[gp], K, dND, filler.Ke.vals);
//			}
//
//			filler.insert();
//		}
//	}
}

struct KSymmFiller
{
	KSymmFiller(HeatTransferKernel &kernel, Kernel::InstanceFiller &filler)
	: stiffness(kernel._stiffness, filler.interval, kernel._stiffness.nodes[filler.interval]),
	  filler(filler) {}

	Kernel::InputParameterIterator stiffness;
	Kernel::InstanceFiller &filler;

	void operator()()
	{
		for (int r = 0; r < stiffness.inc; ++r) {
			for (int c = 0; c < stiffness.inc; ++c, ++stiffness.data) {
				if (r <= c) {
					filler.K[*filler.offset++] += *stiffness.data;
				}
			}
		}
	}

	void operator++() { }
};

struct KFullFiller
{
	KFullFiller(HeatTransferKernel &kernel, Kernel::InstanceFiller &filler)
	: stiffness(kernel._stiffness, filler.interval, kernel._stiffness.nodes[filler.interval]),
	  filler(filler) {}

	Kernel::InputParameterIterator stiffness;
	Kernel::InstanceFiller &filler;

	void operator()()
	{
		for (int r = 0; r < stiffness.inc; ++r) {
			for (int c = 0; c < stiffness.inc; ++c, ++stiffness.data) {
				filler.K[*filler.offset++] += *stiffness.data;
			}
		}
	}

	void operator++() { }
};

void HeatTransferKernel::processElements(const Builder &builder, InstanceFiller &filler)
{
	each_element(KSymmFiller(*this, filler), filler.interval);
}

void HeatTransferKernel::processBoundary(const Builder &builder, size_t rindex, InstanceFiller &filler)
{

}

void HeatTransferKernel::processSolution()
{

}

void HeatTransferKernel::assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, MatrixDense &K, MatrixDense &CD, bool tangentCorrection) const
{

}

struct HeatStiffness
{
	HeatStiffness(HeatTransferKernel &kernel, int interval)
	: w(kernel._w, interval, kernel._w.gps[interval]),
	  det(kernel._detJ, interval, 1),
	  dND(kernel._dND, interval, kernel._dND.dimensions * kernel._dND.nodes[interval]),
	  stiffness(kernel._stiffness, interval, kernel._stiffness.nodes[interval] * kernel._stiffness.nodes[interval])
	{ }

	Kernel::InputParameterIterator w, det, dND;
	Kernel::OuputParameterIterator stiffness;
};

struct Heat2DStiffnessIsotropic: public HeatStiffness
{
	Heat2DStiffnessIsotropic(HeatTransferKernel &kernel, int interval)
	: HeatStiffness(kernel, interval),
	  K(kernel._isotropicK, interval, 1),
	  thickness(kernel._thickness, interval, 1)
	{ }

	Kernel::InputParameterIterator K, thickness;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		KMN2M2N<nodes>(K[0] * thickness[0] * det[0] * w[gpindex], dND.data, stiffness.data);
	}

	void operator++()
	{
		++thickness; ++det; ++dND; ++K;
	}

	void next()
	{
		++stiffness;
	}
};

void HeatTransferKernel::computeStiffness()
{
	_stiffness.forceResize();
	std::fill(_stiffness.data->datatarray().begin(), _stiffness.data->datatarray().end(), 0);

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		each_gp(Heat2DStiffnessIsotropic(*this, i), i);
	}
}

void HeatTransferKernel::computeThickness()
{
	_thickness.smartResize();
	_thickness.eval();
}

struct CoordinatesToNodes
{
	CoordinatesToNodes(HeatTransferKernel &kernel, int interval)
	: procNodes(info::mesh->elements->procNodes->cbegin() + info::mesh->elements->eintervals[interval].begin),
	  ncoordinates(kernel._ncoordinates, interval, kernel._ncoordinates.dimensions)
	{ }

	serializededata<esint, esint>::const_iterator procNodes;
	Kernel::OuputParameterIterator ncoordinates;

	void operator++()
	{
		++procNodes;
	}
};

struct Coordinates2DToNodes: CoordinatesToNodes {
	Coordinates2DToNodes(HeatTransferKernel &kernel, esint interval): CoordinatesToNodes(kernel, interval) {}

	void operator()()
	{
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
			ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
			ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
		}
	}
};

struct Coordinates3DToNodes: CoordinatesToNodes {
	Coordinates3DToNodes(HeatTransferKernel &kernel, esint interval): CoordinatesToNodes(kernel, interval) {}

	void operator()()
	{
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
			ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
			ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
			ncoordinates[2] = info::mesh->nodes->coordinates->datatarray()[*n].z;
		}
	}
};

template<int dimension>
struct ApplyBaseFunctions
{
	ApplyBaseFunctions(HeatTransferKernel &kernel, Kernel::ElementParameterInfo &n, Kernel::ElementParameterInfo &gp, int interval)
	: basis(kernel._N, interval, kernel._N.nodes[interval]),
	  n(n, interval, dimension * n.nodes[interval]),
	  gp(gp, interval, dimension)
	{ }

	Kernel::InputParameterIterator basis, n;
	Kernel::OuputParameterIterator gp;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		NtoGP<nodes, dimension>(basis.data + gpindex * nodes, n.data, gp.data);
	}

	void operator++()
	{
		++gp;
	}

	void next()
	{
		++n;
	}
};


void HeatTransferKernel::computeCoodinates()
{
	_ncoordinates.forceResize();
	_coordinates.forceResize();

	if (info::mesh->dimension == 2) {
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			each_element(Coordinates2DToNodes(*this, i), i);
			each_gp(ApplyBaseFunctions<2>(*this, _ncoordinates, _coordinates, i), i);
		}
	}
	if (info::mesh->dimension == 3) {
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			each_element(Coordinates3DToNodes(*this, i), i);
			each_gp(ApplyBaseFunctions<3>(*this, _ncoordinates, _coordinates, i), i);
		}
	}
}

void HeatTransferKernel::initTemperature()
{
	_ntemperature.forceResize();
	_ninitTemperature.forceResize();
	_temperature.forceResize();
	_initTemperature.forceResize();

	_ninitTemperature.eval();
	_ntemperature.data->operator =(*_ninitTemperature.data);

	move(_ntemperature, temperature);
	initialTemperature->data = temperature->data;
	if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
		// TODO: set
	}
	move(temperature, _ntemperature);

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		each_gp(ApplyBaseFunctions<1>(*this, _ntemperature, _temperature, i), i);
	}
	move(_temperature, _initTemperature);
}

struct Jacobian
{
	Jacobian(HeatTransferKernel &kernel, int interval)
	: coords(kernel._ncoordinates, interval, kernel._ncoordinates.dimensions * kernel._ncoordinates.nodes[interval]),
	  dN(kernel._dN, interval, kernel._dN.dimensions * kernel._dN.nodes[interval]),
	  inv(kernel._invJ, interval, kernel._invJ.dimensions),
	  det(kernel._detJ, interval, 1),
	  dND(kernel._dND, interval, kernel._dND.dimensions * kernel._dND.nodes[interval])
	{ }

	Kernel::InputParameterIterator coords, dN;
	Kernel::OuputParameterIterator inv, det, dND;

	void operator++()
	{
		++inv; ++det; ++dND;
	}

	void next()
	{
		++coords;
	}
};

struct Jacobian2D: public Jacobian
{
	Jacobian2D(HeatTransferKernel &kernel, int interval)
	: Jacobian(kernel, interval) { }

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double jacobian[4] = { 0, 0, 0, 0 };

		for (int n = 0; n < nodes; ++n) {
			jacobian[0] += dN[2 * gpindex * nodes + n + 0 * nodes] * coords[2 * n + 0];
			jacobian[1] += dN[2 * gpindex * nodes + n + 0 * nodes] * coords[2 * n + 1];
			jacobian[2] += dN[2 * gpindex * nodes + n + 1 * nodes] * coords[2 * n + 0];
			jacobian[3] += dN[2 * gpindex * nodes + n + 1 * nodes] * coords[2 * n + 1];
		}

		det[0] = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
		double detJx = 1 / det[0];
		inv[0] =   detJx * jacobian[3];
		inv[1] = - detJx * jacobian[1];
		inv[2] = - detJx * jacobian[2];
		inv[3] =   detJx * jacobian[0];

		M22M2N<nodes>(inv.data, dN.data + 2 * gpindex * nodes, dND.data);
	}
};

struct Jacobian3D: public Jacobian
{
	Jacobian3D(HeatTransferKernel &kernel, int interval)
	: Jacobian(kernel, interval) { }

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		double jacobian[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

		for (int n = 0; n < nodes; ++n) {
			jacobian[0] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 0];
			jacobian[1] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 1];
			jacobian[2] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 2];
			jacobian[3] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 0];
			jacobian[4] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 1];
			jacobian[5] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 2];
			jacobian[6] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 0];
			jacobian[7] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 1];
			jacobian[8] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 2];
		}
		det[0] =
				+ jacobian[0] * jacobian[4] * jacobian[8]
				+ jacobian[1] * jacobian[5] * jacobian[6]
				+ jacobian[2] * jacobian[3] * jacobian[7]
				- jacobian[2] * jacobian[4] * jacobian[6]
				- jacobian[1] * jacobian[3] * jacobian[8]
				- jacobian[0] * jacobian[5] * jacobian[7];

		double detJx = 1 / det[0];
		inv[0] = detJx * ( jacobian[8] * jacobian[4] - jacobian[7] * jacobian[5]);
		inv[1] = detJx * (-jacobian[8] * jacobian[1] + jacobian[7] * jacobian[2]);
		inv[2] = detJx * ( jacobian[5] * jacobian[1] - jacobian[4] * jacobian[2]);
		inv[3] = detJx * (-jacobian[8] * jacobian[3] + jacobian[6] * jacobian[5]);
		inv[4] = detJx * ( jacobian[8] * jacobian[0] - jacobian[6] * jacobian[2]);
		inv[5] = detJx * (-jacobian[5] * jacobian[0] + jacobian[3] * jacobian[2]);
		inv[6] = detJx * ( jacobian[7] * jacobian[3] - jacobian[6] * jacobian[4]);
		inv[7] = detJx * (-jacobian[7] * jacobian[0] + jacobian[6] * jacobian[1]);
		inv[8] = detJx * ( jacobian[4] * jacobian[0] - jacobian[3] * jacobian[1]);

		M33M3N<nodes>(inv.data, dN.data + 3 * gpindex * nodes, dND.data);
	}
};


void HeatTransferKernel::computeJacobian()
{
	_invJ.forceResize();
	_detJ.forceResize();
	_dND.forceResize();

	if (info::mesh->dimension == 2) {
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			each_gp(Jacobian2D(*this, i), i);
		}
	}
	if (info::mesh->dimension == 3) {
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			each_gp(Jacobian3D(*this, i), i);
		}
	}
}

struct Advection
{
//	Advection(HeatTransferKernel &kernel, esint interval)
//	: dND(kernel._dND, interval),
//	  motion(kernel._translationMotion, interval),
//	  dens(kernel._dens, interval),
//	  hc(kernel._hc, interval),
//	  advection(kernel._advection, interval),
//	  u(kernel._u, interval)
//	{ }
//
//	Kernel::InputParameterIterator dND, motion, dens, hc;
//	Kernel::OuputParameterIterator advection, u;
//
//	template<int dimension, int N, int GP>
//	void operator()()
//	{
//		for (int d = 0; d < dimension; ++d) {
//			u[d] = motion[d] * dens[0] * hc[0];
//		}
//		M12M2N<N>(u.data, dND.data, advection.data);
//	}
//
//	void operator++()
//	{
//		++dND; ++motion; ++dens, ++hc;
//		++advection, ++u;
//	}
};

void HeatTransferKernel::computeAdvection()
{
//	if (_translationMotion.constness == ParameterInfo::Status::EMPTY) {
//		return;
//	}
//
//	resize(_dens, [&] (const ElementsInterval &ei) -> int { return Mesh::edata[ei.code].N->size(); });
//	resize(_hc, [&] (const ElementsInterval &ei) -> int { return Mesh::edata[ei.code].N->size(); });
//	resize(_translationMotion, [&] (const ElementsInterval &ei) -> int { return info::mesh->dimension * Mesh::edata[ei.code].N->size(); });
//	resize(_u, [&] (const ElementsInterval &ei) -> int { return info::mesh->dimension * Mesh::edata[ei.code].N->size(); });
//	resize(_advection, [&] (const ElementsInterval &ei) -> int { return Mesh::edata[ei.code].nodes * Mesh::edata[ei.code].N->size(); });
//
//	evaluate(_dens);
//	evaluate(_hc);
//	evaluate(_translationMotion);

//	each_gp<Advection>(*this);
}

struct ConductivityIsotropic {

	ConductivityIsotropic(HeatTransferKernel &kernel, int interval)
	: conductivity(kernel._conductivity[0], interval, 1),
	  K(kernel._isotropicK, interval, 1)
	{ }

	Kernel::InputParameterIterator conductivity;
	Kernel::OuputParameterIterator K;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		*K.data = *conductivity.data;
	}

	void operator++() { ++K; ++conductivity; }
	void next() {}
};

void HeatTransferKernel::computeConductivity()
{
	_csCartesian.smartResize();
	_csSpherical.smartResize();
	_csCylindrical.smartResize();

	_conductivity[0].smartResize();
	_conductivity[1].smartResize();
	_conductivity[2].smartResize();
	_conductivity[3].smartResize();

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		int material = info::mesh->elements->eintervals[i].material;
		switch (info::mesh->materials[material]->thermal_conductivity.model) {
		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: _K.isconstant[i] &= _conductivity[3].constantInterval(i); break;
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:   _K.isconstant[i] &= _conductivity[2].constantInterval(i); break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:    _K.isconstant[i] &= _conductivity[1].constantInterval(i); break;
		case ThermalConductivityConfiguration::MODEL::ISOTROPIC: _isotropicK.isconstant[i] &= _conductivity[0].constantInterval(i); break;
		}
	}

	_isotropicK.smartResize();
	_K.smartResize();

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		int material = info::mesh->elements->eintervals[i].material;
		switch (info::mesh->materials[material]->thermal_conductivity.model) {
		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
			break;
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:
			break;
		case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
			_conductivity[0].eval();
			each_gp(ConductivityIsotropic(*this, i), i);
			break;
		}
	}


//	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
//		int material = info::mesh->elements->eintervals[i].material;
//		if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//			_K.isconstant[i] = _csCartesian.constantInterval(i);
//		}
//		if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL) {
//			_K.isconstant[i] = false;
//			_K.addGeneralInput(_coordinates);
//		}
//		if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CYLINDRICAL) {
//			_K.isconstant[i] = false;
//			_K.addGeneralInput(_coordinates);
//		}
//
//		switch (info::mesh->materials[material]->thermal_conductivity.model) {
//		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: _K.isconstant[i] &= _conductivity[3].constantInterval(i); break;
//		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:   _K.isconstant[i] &= _conductivity[2].constantInterval(i); break;
//		case ThermalConductivityConfiguration::MODEL::DIAGONAL:    _K.isconstant[i] &= _conductivity[1].constantInterval(i); break;
//		case ThermalConductivityConfiguration::MODEL::ISOTROPIC: _isotropicK.isconstant[i] &= _conductivity[0].constantInterval(i); break;
//		}
//	}
//
//	resize(_K, [&] (const ElementsInterval &ei) -> int { return info::mesh->dimension * info::mesh->dimension * Mesh::edata[ei.code].N->size(); });
//	resize(_isotropicK, [&] (const ElementsInterval &ei) -> int { return Mesh::edata[ei.code].N->size(); });
//
//	evaluate(_csCartesian);
//	evaluate(_csSpherical);
//	evaluate(_csCylindrical);
//	evaluate(_conductivity[0]);
//	evaluate(_conductivity[1]);
//	evaluate(_conductivity[2]);
//	evaluate(_conductivity[3]);
//
//	if (info::mesh->dimension == 2) {
//		computeConductivity2D();
//	}
//	if (info::mesh->dimension == 3) {
//		computeConductivity3D();
//	}
}

void HeatTransferKernel::computeConductivity2D()
{
//	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
//		int material = info::mesh->elements->eintervals[i].material;
//		double *c = NULL; int cinc = 0;
//		auto begin = (_K.data->begin() + i)->begin();
//		auto end = (_K.data->begin() + i)->end();
//
//		switch (info::mesh->materials[material]->thermal_conductivity.model) {
//		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: {
//			c = (_conductivity[3].data->begin() + i)->begin(); cinc = _conductivity[3].constantInterval(i) ? 0 : 4;
//		} break;
//		case ThermalConductivityConfiguration::MODEL::SYMMETRIC: {
//			c = (_conductivity[2].data->begin() + i)->begin(); cinc = _conductivity[2].constantInterval(i) ? 0 : 3;
//		} break;
//		case ThermalConductivityConfiguration::MODEL::DIAGONAL: {
//			c = (_conductivity[1].data->begin() + i)->begin(); cinc = _conductivity[1].constantInterval(i) ? 0 : 2;
//		} break;
//		case ThermalConductivityConfiguration::MODEL::ISOTROPIC: {
//			c = (_conductivity[0].data->begin() + i)->begin(); cinc = _conductivity[0].constantInterval(i) ? 0 : 1;
//			begin = (_isotropicK.data->begin() + i)->begin();
//			end = (_isotropicK.data->begin() + i)->end();
//		} break;
//		}
//
//		if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//			double value = 0;
//			if (_csCartesian.isconstant[i]) {
//				value = _csCartesian.ivalues[i];
//			} else {
//				Evaluator::Params params;
//				_csCartesian.insertParameterValues(params, i, 0);
//				value = _csCartesian.evaluator[i]->eval(params);
//			}
//			double cos = std::cos(M_PI * value / 180);
//			double sin = std::sin(M_PI * value / 180);
//
//			switch (info::mesh->materials[material]->thermal_conductivity.model) {
//			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: {
//				for (auto k = begin; k != end; k += 4, c += cinc) {
//					rotate2Dfull(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: {
//				for (auto k = begin; k != end; k += 4, c += cinc) {
//					rotate2Dsym(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::DIAGONAL: {
//				for (auto k = begin; k != end; k += 4, c += cinc) {
//					rotate2Ddiag(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: {
//				for (auto k = begin; k != end; ++k, c += cinc) {
//					k[0] = c[0];
//				}
//			} break;
//			}
//		}
//		if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CYLINDRICAL) {
//			std::vector<double, initless_allocator<double> > centerx, centery;
//			bool isotropic = info::mesh->materials[material]->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
//			if (_csCylindrical.isconstant[i + 1] || isotropic) {
//				centerx.push_back(_csCylindrical.ivalues[2 * i + 0]);
//			} else {
//				centerx.resize((_coordinates.data->begin() + i)->size() / 2);
//				Evaluator::Params params;
//				_csCylindrical.insertParameterValues(params, i, 1);
//				_csCylindrical.evaluator[2 *i + 0]->evalVector(centerx.size(), 1, params, centerx.data());
//			}
//			if (_csCylindrical.isconstant[2 *i + 1] || isotropic) {
//				centery.push_back(_csCylindrical.ivalues[2 * i + 1]);
//			} else {
//				centery.resize((_coordinates.data->begin() + i)->size() / 2);
//				Evaluator::Params params;
//				_csCylindrical.insertParameterValues(params, i, 1);
//				_csCylindrical.evaluator[2 * i + 1]->evalVector(centery.size(), 1, params, centery.data());
//			}
//			double *ox = centerx.data(), *oy = centery.data();
//			double *coo = (_coordinates.data->begin() + i)->begin();
//			int oxinc = _csCylindrical.isconstant[2 * i + 0] ? 0 : 1, oyinc = _csCylindrical.isconstant[2 * i + 1] ? 0 : 1;
//
//			auto t = [] (double &cos, double &sin, const double *coo, const double *ox, const double *oy) {
//				double rotation = std::atan2((coo[1] - *oy), (coo[0] - *ox));
//				cos = std::cos(rotation); sin = std::sin(rotation);
//			};
//
//			switch (info::mesh->materials[material]->thermal_conductivity.model) {
//			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: {
//				for (auto k = begin; k != end; k += 4, c += cinc, coo += 2, ox += oxinc, oy += oyinc) {
//					double cos, sin; t(cos, sin, coo, ox, oy); rotate2Dfull(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: {
//				for (auto k = begin; k != end; k += 4, c += cinc, coo += 2, ox += oxinc, oy += oyinc) {
//					double cos, sin; t(cos, sin, coo, ox, oy); rotate2Dsym(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::DIAGONAL: {
//				for (auto k = begin; k != end; k += 4, c += cinc, coo += 2, ox += oxinc, oy += oyinc) {
//					double cos, sin; t(cos, sin, coo, ox, oy); rotate2Ddiag(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: {
//				for (auto k = begin; k != end; ++k, c += cinc) {
//					k[0] = c[0];
//				}
//			} break;
//			}
//		}
//	}
}

void HeatTransferKernel::computeConductivity3D()
{
//	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
//		int material = info::mesh->elements->eintervals[i].material;
//		double *c = NULL, cos[3], sin[3]; int cinc = 0;
//		auto begin = (_K.data->begin() + i)->begin();
//		auto end = (_K.data->begin() + i)->end();
//
//		switch (info::mesh->materials[material]->thermal_conductivity.model) {
//		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: {
//			c = (_conductivity[3].data->begin() + i)->begin(); cinc = _conductivity[3].constantInterval(i) ? 0 : 9;
//		} break;
//		case ThermalConductivityConfiguration::MODEL::SYMMETRIC: {
//			c = (_conductivity[2].data->begin() + i)->begin(); cinc = _conductivity[2].constantInterval(i) ? 0 : 6;
//		} break;
//		case ThermalConductivityConfiguration::MODEL::DIAGONAL: {
//			c = (_conductivity[1].data->begin() + i)->begin(); cinc = _conductivity[1].constantInterval(i) ? 0 : 3;
//		} break;
//		case ThermalConductivityConfiguration::MODEL::ISOTROPIC: {
//			c = (_conductivity[0].data->begin() + i)->begin(); cinc = _conductivity[0].constantInterval(i) ? 0 : 1;
//			begin = (_isotropicK.data->begin() + i)->begin();
//			end = (_isotropicK.data->begin() + i)->end();
//		} break;
//		}
//
//		if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//			auto eval = [&] (int d) {
//				if (_csCartesian.isconstant[3 * i + d]) {
//					cos[d] = std::cos(M_PI * _csCartesian.ivalues[i] / 180);
//					sin[d] = std::sin(M_PI * _csCartesian.ivalues[i] / 180);
//				} else {
//					Evaluator::Params params;
//					_csCartesian.insertParameterValues(params, i, d);
//					double value = _csCartesian.evaluator[i]->eval(params);
//					cos[d] = std::cos(M_PI * value / 180);
//					sin[d] = std::sin(M_PI * value / 180);
//				}
//			};
//			eval(0); eval(1); eval(2);
//
//			switch (info::mesh->materials[material]->thermal_conductivity.model) {
//			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: {
//				for (auto k = begin; k != end; k += 9, c += cinc) {
//					rotate3Dfull(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: {
//				for (auto k = begin; k != end; k += 9, c += cinc) {
//					rotate3Dsym(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::DIAGONAL: {
//				for (auto k = begin; k != end; k += 9, c += cinc) {
//					rotate3Ddiag(cos, sin, c, k);
//				}
//			} break;
//			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: {
//				for (auto k = begin; k != end; ++k, c += cinc) {
//					k[0] = c[0];
//				}
//			} break;
//			}
//		} else {
//			std::vector<double, initless_allocator<double> > centerx, centery, centerz;
//			auto eval = [&] (ParameterInfo &info, std::vector<double, initless_allocator<double> > &values, int d) {
//				bool isotropic = info::mesh->materials[material]->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC;
//				if (info.isconstant[i + 1] || isotropic) {
//					values.push_back(info.ivalues[2 * i + d]);
//				} else {
//					values.resize((_coordinates.data->begin() + i)->size() / 2);
//					Evaluator::Params params;
//					info.insertParameterValues(params, i, d);
//					info.evaluator[2 * i + 0]->evalVector(values.size(), 1, params, values.data());
//				}
//			};
//
//			if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::SPHERICAL) {
//				eval(_csSpherical, centerx, 0);
//				eval(_csSpherical, centery, 1);
//				eval(_csSpherical, centerz, 2);
//				auto begin = (_K.data->begin() + i)->begin();
//				auto end = (_K.data->begin() + i)->end();
//				double *ox = centerx.data();
//				double *oy = centery.data();
//				double *oz = centerz.data();
//				double *coo = (_coordinates.data->begin() + i)->begin();
//				int oxinc = _csSpherical.isconstant[3 * i + 0] ? 0 : 1;
//				int oyinc = _csSpherical.isconstant[3 * i + 1] ? 0 : 1;
//				int ozinc = _csSpherical.isconstant[3 * i + 2] ? 0 : 1;
//
//				cos[0] = 1.0;
//				sin[0] = 0.0;
//				auto t = [] (double *cos, double *sin, const double *coo, const double *ox, const double *oy, const double *oz) {
//					double azimut = std::atan2((coo[1] - *oy), (coo[0] - *ox));
//					double elevation = 0.0;
//					if ((coo[0] - *ox) * (coo[0] - *ox) + (coo[1] - *oy) * (coo[1] - *oy) + (coo[2] - *oz) * (coo[2] - *oz) >= 1e-15) {
//						elevation = std::atan2(std::sqrt((coo[2] - *oz) * (coo[2] - *oz) + (coo[0] - *ox) * (coo[0] - *ox)), (coo[1] - *oy));
//					}
//					cos[1] = std::cos(elevation);
//					cos[2] = std::cos(azimut);
//					sin[1] = std::sin(elevation);
//					sin[2] = std::sin(azimut);
//				};
//
//				switch (info::mesh->materials[material]->thermal_conductivity.model) {
//				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: {
//					for (auto k = begin; k != end; k += 9, c += cinc, coo += 3, ox += oxinc, oy += oyinc, oz += ozinc) {
//						t(cos, sin, coo, ox, oy, oz); rotate3Dfull(cos, sin, c, k);
//					}
//				} break;
//				case ThermalConductivityConfiguration::MODEL::SYMMETRIC: {
//					for (auto k = begin; k != end; k += 9, c += cinc, coo += 3, ox += oxinc, oy += oyinc, oz += ozinc) {
//						t(cos, sin, coo, ox, oy, oz); rotate3Dsym(cos, sin, c, k);
//					}
//				} break;
//				case ThermalConductivityConfiguration::MODEL::DIAGONAL: {
//					for (auto k = begin; k != end; k += 9, c += cinc, coo += 3, ox += oxinc, oy += oyinc, oz += ozinc) {
//						t(cos, sin, coo, ox, oy, oz); rotate3Ddiag(cos, sin, c, k);
//					}
//				} break;
//				case ThermalConductivityConfiguration::MODEL::ISOTROPIC: {
//					for (auto k = begin; k != end; ++k, c += cinc) {
//						k[0] = c[0];
//					}
//				} break;
//				}
//			}
//
//			if (info::mesh->materials[material]->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CYLINDRICAL) {
//				eval(_csCylindrical, centerx, 0);
//				eval(_csCylindrical, centery, 1);
//				auto begin = (_K.data->begin() + i)->begin();
//				auto end = (_K.data->begin() + i)->end();
//				double *ox = centerx.data();
//				double *oy = centery.data();
//				double *coo = (_coordinates.data->begin() + i)->begin();
//				int oxinc = _csCylindrical.isconstant[2 * i + 0] ? 0 : 1;
//				int oyinc = _csCylindrical.isconstant[2 * i + 1] ? 0 : 1;
//
//				cos[0] = 1.0;
//				cos[1] = 1.0;
//				sin[0] = 0.0;
//				sin[1] = 0.0;
//				auto t = [] (double *cos, double *sin, const double *coo, const double *ox, const double *oy) {
//					double rotation = std::atan2((coo[1] - *oy), (coo[0] - *ox));
//					cos[2] = std::cos(rotation);
//					sin[2] = std::sin(rotation);
//				};
//
//				switch (info::mesh->materials[material]->thermal_conductivity.model) {
//				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: {
//					for (auto k = begin; k != end; k += 9, c += cinc, coo += 3, ox += oxinc, oy += oyinc) {
//						t(cos, sin, coo, ox, oy); rotate3Dfull(cos, sin, c, k);
//					}
//				} break;
//				case ThermalConductivityConfiguration::MODEL::SYMMETRIC: {
//					for (auto k = begin; k != end; k += 9, c += cinc, coo += 3, ox += oxinc, oy += oyinc) {
//						t(cos, sin, coo, ox, oy); rotate3Dsym(cos, sin, c, k);
//					}
//				} break;
//				case ThermalConductivityConfiguration::MODEL::DIAGONAL: {
//					for (auto k = begin; k != end; k += 9, c += cinc, coo += 3, ox += oxinc, oy += oyinc) {
//						t(cos, sin, coo, ox, oy); rotate3Ddiag(cos, sin, c, k);
//					}
//				} break;
//				case ThermalConductivityConfiguration::MODEL::ISOTROPIC: {
//					for (auto k = begin; k != end; ++k, c += cinc) {
//						k[0] = c[0];
//					}
//				} break;
//				}
//			}
//		}
//	}
}











