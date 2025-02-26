
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "math/simd/simd.h"

namespace espreso {

template <class NGP, class Operator>
static inline void iterate_elements_gps(Operator op)
{
	if (op.update) {
		if (op.isconst) {
			switch (info::mesh->elements->eintervals[op.interval].code) {
			case static_cast<size_t>(Element::CODE::LINE2):
				for (size_t gp = 0; gp < NGP::LINE2; ++gp) {
					op.template operator()<2, NGP::LINE2>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::LINE3):
				for (size_t gp = 0; gp < NGP::LINE3; ++gp) {
					op.template operator()<3, NGP::LINE3>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE3):
				for (size_t gp = 0; gp < NGP::TRIANGLE3; ++gp) {
					op.template operator()<3, NGP::TRIANGLE3>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE6):
				for (size_t gp = 0; gp < NGP::TRIANGLE6; ++gp) {
					op.template operator()<6, NGP::TRIANGLE6>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE4):
				for (size_t gp = 0; gp < NGP::SQUARE4; ++gp) {
					op.template operator()<4, NGP::SQUARE4>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE8):
				for (size_t gp = 0; gp < NGP::SQUARE8; ++gp) {
					op.template operator()<8, NGP::SQUARE8>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA4):
				for (size_t gp = 0; gp < NGP::TETRA4; ++gp) {
					op.template operator()<4, NGP::TETRA4>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA10):
				for (size_t gp = 0; gp < NGP::TETRA10; ++gp) {
					op.template operator()<10, NGP::TETRA10>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID5):
				for (size_t gp = 0; gp < NGP::PYRAMID5; ++gp) {
					op.template operator()<5, NGP::PYRAMID5>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID13):
				for (size_t gp = 0; gp < NGP::PYRAMID13; ++gp) {
					op.template operator()<13, NGP::PYRAMID13>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA6):
				for (size_t gp = 0; gp < NGP::PRISMA6; ++gp) {
					op.template operator()<6, NGP::PRISMA6>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA15):
				for (size_t gp = 0; gp < NGP::PRISMA15; ++gp) {
					op.template operator()<15, NGP::PRISMA15>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA8):
				for (size_t gp = 0; gp < NGP::HEXA8; ++gp) {
					op.template operator()<8, NGP::HEXA8>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA20):
				for (size_t gp = 0; gp < NGP::HEXA20; ++gp) {
					op.template operator()<20, NGP::HEXA20>(gp);
				}
			break;
			default: break;
			}
		} else {
			switch (info::mesh->elements->eintervals[op.interval].code) {
			case static_cast<size_t>(Element::CODE::LINE2):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::LINE2; ++gp) {
						op.template operator()<2, NGP::LINE2>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::LINE3):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::LINE3; ++gp) {
						op.template operator()<3, NGP::LINE3>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE3):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TRIANGLE3; ++gp) {
						op.template operator()<3, NGP::TRIANGLE3>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE6):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TRIANGLE6; ++gp) {
						op.template operator()<6, NGP::TRIANGLE6>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE4):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::SQUARE4; ++gp) {
						op.template operator()<4, NGP::SQUARE4>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE8):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::SQUARE8; ++gp) {
						op.template operator()<8, NGP::SQUARE8>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA4):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TETRA4; ++gp) {
						op.template operator()<4, NGP::TETRA4>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA10):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TETRA10; ++gp) {
						op.template operator()<10, NGP::TETRA10>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID5):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PYRAMID5; ++gp) {
						op.template operator()<5, NGP::PYRAMID5>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID13):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PYRAMID13; ++gp) {
						op.template operator()<13, NGP::PYRAMID13>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA6):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PRISMA6; ++gp) {
						op.template operator()<6, NGP::PRISMA6>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA15):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PRISMA15; ++gp) {
						op.template operator()<15, NGP::PRISMA15>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA8):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::HEXA8; ++gp) {
						op.template operator()<8, NGP::HEXA8>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA20):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::HEXA20; ++gp) {
						op.template operator()<20, NGP::HEXA20>(gp);
					}
				}
			break;
			default: break;
			}
		}
	}
}

template <class NGP, class Operator>
static inline void iterate_elements_gps_simd(Operator op)
{
	if (op.update) {
		if (op.isconst) {
			switch (info::mesh->elements->eintervals[op.interval].code) {
			case static_cast<size_t>(Element::CODE::LINE2):
				for (size_t gp = 0; gp < NGP::LINE2; ++gp) {
					op.template operator()<2, NGP::LINE2>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::LINE3):
				for (size_t gp = 0; gp < NGP::LINE3; ++gp) {
					op.template operator()<3, NGP::LINE3>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE3):
				for (size_t gp = 0; gp < NGP::TRIANGLE3; ++gp) {
					op.template operator()<3, NGP::TRIANGLE3>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE6):
				for (size_t gp = 0; gp < NGP::TRIANGLE6; ++gp) {
					op.template operator()<6, NGP::TRIANGLE6>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE4):
				for (size_t gp = 0; gp < NGP::SQUARE4; ++gp) {
					op.template operator()<4, NGP::SQUARE4>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE8):
				for (size_t gp = 0; gp < NGP::SQUARE8; ++gp) {
					op.template operator()<8, NGP::SQUARE8>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA4):
				for (size_t gp = 0; gp < NGP::TETRA4; ++gp) {
					op.template operator()<4, NGP::TETRA4>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA10):
				for (size_t gp = 0; gp < NGP::TETRA10; ++gp) {
					op.template operator()<10, NGP::TETRA10>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID5):
				for (size_t gp = 0; gp < NGP::PYRAMID5; ++gp) {
					op.template operator()<5, NGP::PYRAMID5>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID13):
				for (size_t gp = 0; gp < NGP::PYRAMID13; ++gp) {
					op.template operator()<13, NGP::PYRAMID13>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA6):
				for (size_t gp = 0; gp < NGP::PRISMA6; ++gp) {
					op.template operator()<6, NGP::PRISMA6>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA15):
				for (size_t gp = 0; gp < NGP::PRISMA15; ++gp) {
					op.template operator()<15, NGP::PRISMA15>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA8):
				for (size_t gp = 0; gp < NGP::HEXA8; ++gp) {
					op.template operator()<8, NGP::HEXA8>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA20):
				for (size_t gp = 0; gp < NGP::HEXA20; ++gp) {
					op.template operator()<20, NGP::HEXA20>(gp);
				}
			break;
			default: break;
			}
		} else {
			switch (info::mesh->elements->eintervals[op.interval].code) {
			case static_cast<size_t>(Element::CODE::LINE2):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::LINE2; ++gp) {
						op.template operator()<2, NGP::LINE2>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::LINE3):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::LINE3; ++gp) {
						op.template operator()<3, NGP::LINE3>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE3):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::TRIANGLE3; ++gp) {
						op.template operator()<3, NGP::TRIANGLE3>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE6):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::TRIANGLE6; ++gp) {
						op.template operator()<6, NGP::TRIANGLE6>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE4):
				__asm volatile("# Loop calling Stiffness2DHeatSimd::operator() method start");
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {	
					for (size_t gp = 0; gp < NGP::SQUARE4; ++gp) {
						op.template operator()<4, NGP::SQUARE4>(gp);
					}
			}
				__asm volatile("#Loop calling  Stiffness2DHeatSimd::operator() method end");
			break;
			case static_cast<size_t>(Element::CODE::SQUARE8):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::SQUARE8; ++gp) {
						op.template operator()<8, NGP::SQUARE8>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA4):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::TETRA4; ++gp) {
						op.template operator()<4, NGP::TETRA4>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA10):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::TETRA10; ++gp) {
						op.template operator()<10, NGP::TETRA10>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID5):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end;i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::PYRAMID5; ++gp) {
						op.template operator()<5, NGP::PYRAMID5>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID13):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::PYRAMID13; ++gp) {
						op.template operator()<13, NGP::PYRAMID13>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA6):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::PRISMA6; ++gp) {
						op.template operator()<6, NGP::PRISMA6>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA15):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::PRISMA15; ++gp) {
						op.template operator()<15, NGP::PRISMA15>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA8):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::HEXA8; ++gp) {
						op.template operator()<8, NGP::HEXA8>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA20):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
					for (size_t gp = 0; gp < NGP::HEXA20; ++gp) {
						op.template operator()<20, NGP::HEXA20>(gp);
					}
				}
			break;
			default: break;
			}
		}
	}
}

template <class Operator>
static inline void iterate_elements(Operator op)
{
	if (op.update) {
		if (op.isconst) {
			op();
		} else {
			for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
				op();
			}
		}
	}
}

template <class Operator>
static inline void iterate_elements_simd(Operator op)
{
	if (op.update) {
		if (op.isconst) {
			op();
		} else {
			for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; i+=SIMD::size, op+=SIMD::size) {
				op();
			}
		}
	}
}

template <class NGP, class Operator>
static inline void iterate_boundary_gps(Operator op, size_t region)
{
	if (op.update) {
		if (op.isconst) {
			switch (info::mesh->boundaryRegions[region]->eintervals[op.interval].code) {
			case static_cast<size_t>(Element::CODE::LINE2):
				for (size_t gp = 0; gp < NGP::LINE2; ++gp) {
					op.template operator()<2, NGP::LINE2>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::LINE3):
			for (size_t gp = 0; gp < NGP::LINE3; ++gp) {
				op.template operator()<3, NGP::LINE3>(gp);
			}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE3):
				for (size_t gp = 0; gp < NGP::TRIANGLE3; ++gp) {
					op.template operator()<3, NGP::TRIANGLE3>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE6):
				for (size_t gp = 0; gp < NGP::TRIANGLE6; ++gp) {
					op.template operator()<6, NGP::TRIANGLE6>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE4):
				for (size_t gp = 0; gp < NGP::SQUARE4; ++gp) {
					op.template operator()<4, NGP::SQUARE4>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE8):
				for (size_t gp = 0; gp < NGP::SQUARE8; ++gp) {
					op.template operator()<8, NGP::SQUARE8>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA4):
				for (size_t gp = 0; gp < NGP::TETRA4; ++gp) {
					op.template operator()<4, NGP::TETRA4>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA10):
				for (size_t gp = 0; gp < NGP::TETRA10; ++gp) {
					op.template operator()<10, NGP::TETRA10>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID5):
				for (size_t gp = 0; gp < NGP::PYRAMID5; ++gp) {
					op.template operator()<5, NGP::PYRAMID5>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID13):
				for (size_t gp = 0; gp < NGP::PYRAMID13; ++gp) {
					op.template operator()<13, NGP::PYRAMID13>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA6):
				for (size_t gp = 0; gp < NGP::PRISMA6; ++gp) {
					op.template operator()<6, NGP::PRISMA6>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA15):
				for (size_t gp = 0; gp < NGP::PRISMA15; ++gp) {
					op.template operator()<15, NGP::PRISMA15>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA8):
				for (size_t gp = 0; gp < NGP::HEXA8; ++gp) {
					op.template operator()<8, NGP::HEXA8>(gp);
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA20):
				for (size_t gp = 0; gp < NGP::HEXA20; ++gp) {
					op.template operator()<20, NGP::HEXA20>(gp);
				}
			break;
			default: break;
			}
		} else {
			switch (info::mesh->boundaryRegions[region]->eintervals[op.interval].code) {
			case static_cast<size_t>(Element::CODE::LINE2):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::LINE2; ++gp) {
						op.template operator()<2, NGP::LINE2>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::LINE3):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::LINE3; ++gp) {
						op.template operator()<3, NGP::LINE3>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE3):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TRIANGLE3; ++gp) {
						op.template operator()<3, NGP::TRIANGLE3>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TRIANGLE6):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TRIANGLE6; ++gp) {
						op.template operator()<6, NGP::TRIANGLE6>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE4):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::SQUARE4; ++gp) {
						op.template operator()<4, NGP::SQUARE4>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::SQUARE8):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::SQUARE8; ++gp) {
						op.template operator()<8, NGP::SQUARE8>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA4):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TETRA4; ++gp) {
						op.template operator()<4, NGP::TETRA4>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::TETRA10):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::TETRA10; ++gp) {
						op.template operator()<10, NGP::TETRA10>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID5):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PYRAMID5; ++gp) {
						op.template operator()<5, NGP::PYRAMID5>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PYRAMID13):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PYRAMID13; ++gp) {
						op.template operator()<13, NGP::PYRAMID13>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA6):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PRISMA6; ++gp) {
						op.template operator()<6, NGP::PRISMA6>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::PRISMA15):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::PRISMA15; ++gp) {
						op.template operator()<15, NGP::PRISMA15>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA8):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::HEXA8; ++gp) {
						op.template operator()<8, NGP::HEXA8>(gp);
					}
				}
			break;
			case static_cast<size_t>(Element::CODE::HEXA20):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (size_t gp = 0; gp < NGP::HEXA20; ++gp) {
						op.template operator()<20, NGP::HEXA20>(gp);
					}
				}
			break;
			default: break;
			}
		}
	}
}

template <class Operator>
static inline void iterate_boundary(Operator op, size_t region)
{
	if (op.update) {
		if (op.isconst) {
			op();
		} else {
			for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
				op();
			}
		}
	}
}


}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_ */
