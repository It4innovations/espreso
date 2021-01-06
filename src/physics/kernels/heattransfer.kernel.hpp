
#ifndef SRC_PHYSICS_KERNELS_HEATTRANSFER_KERNEL_HPP_
#define SRC_PHYSICS_KERNELS_HEATTRANSFER_KERNEL_HPP_

#include "heattransfer.kernel.h"
#include "math.kernel.hpp"

namespace espreso {

template <class Iterator>
static inline void each_gp(Iterator it, int interval)
{
	switch (info::mesh->elements->eintervals[interval].code) {
	case static_cast<int>(Element::CODE::TRIANGLE3):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_TRIANGLE3; ++gp, ++it) {
				it.template operator()<3, HeatTransferKernel::GP_TRIANGLE3>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::TRIANGLE6):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_TRIANGLE6; ++gp, ++it) {
				it.template operator()<6, HeatTransferKernel::GP_TRIANGLE6>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::SQUARE4):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_SQUARE4; ++gp, ++it) {
				it.template operator()<4, HeatTransferKernel::GP_SQUARE4>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::SQUARE8):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_SQUARE8; ++gp, ++it) {
				it.template operator()<8, HeatTransferKernel::GP_SQUARE8>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::TETRA4):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_TETRA4; ++gp, ++it) {
				it.template operator()<4, HeatTransferKernel::GP_TETRA4>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::TETRA10):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_TETRA10; ++gp, ++it) {
				it.template operator()<10, HeatTransferKernel::GP_TETRA10>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::PYRAMID5):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_PYRAMID5; ++gp, ++it) {
				it.template operator()<5, HeatTransferKernel::GP_PYRAMID5>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::PYRAMID13):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_PYRAMID13; ++gp, ++it) {
				it.template operator()<13, HeatTransferKernel::GP_PYRAMID13>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::PRISMA6):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_PRISMA6; ++gp, ++it) {
				it.template operator()<6, HeatTransferKernel::GP_PRISMA6>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::PRISMA15):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_PRISMA15; ++gp, ++it) {
				it.template operator()<15, HeatTransferKernel::GP_PRISMA15>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::HEXA8):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_HEXA8; ++gp, ++it) {
				it.template operator()<8, HeatTransferKernel::GP_HEXA8>(gp);
			}
		}
	break;
	case static_cast<int>(Element::CODE::HEXA20):
		for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, it.next()) {
			for (int gp = 0; gp < HeatTransferKernel::GP_HEXA20; ++gp, ++it) {
				it.template operator()<20, HeatTransferKernel::GP_HEXA20>(gp);
			}
		}
	break;
	default: break;
	}
}

template <class Iterator, class Kernel>
static inline void each_gp(Kernel &kernel)
{
	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		Iterator it(kernel, i);
		each_gp(it, i);
	}
}

template <class Iterator>
static inline void each_element(Iterator it, int interval)
{
	for (esint i = info::mesh->elements->eintervals[interval].begin; i < info::mesh->elements->eintervals[interval].end; ++i, ++it) {
		it();
	}
}

template <class Iterator, class Kernel>
static inline void each_element(Kernel &kernel)
{
	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		each_element(Iterator(kernel, i));
	}
}

}

#endif /* SRC_PHYSICS_KERNELS_HEATTRANSFER_KERNEL_HPP_ */
