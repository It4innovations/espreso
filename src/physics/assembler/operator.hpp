
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template <class NGP, class Operator>
static inline void iterate_elements_gps(Operator op)
{
	if (op.update) {
		if (op.isconst) {
			if (Operator::print) printf("\tOP::ELEMENT::%d::GP::CONSTANT::%s\n", op.interval, op.name());
			switch (info::mesh->elements->eintervals[op.interval].code) {
			case static_cast<int>(Element::CODE::LINE2):
				for (int gp = 0; gp < NGP::LINE2; ++gp) {
					op.template operator()<2, NGP::LINE2>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::LINE3):
				for (int gp = 0; gp < NGP::LINE3; ++gp) {
					op.template operator()<3, NGP::LINE3>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE3):
				for (int gp = 0; gp < NGP::TRIANGLE3; ++gp) {
					op.template operator()<3, NGP::TRIANGLE3>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE6):
				for (int gp = 0; gp < NGP::TRIANGLE6; ++gp) {
					op.template operator()<6, NGP::TRIANGLE6>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE4):
				for (int gp = 0; gp < NGP::SQUARE4; ++gp) {
					op.template operator()<4, NGP::SQUARE4>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE8):
				for (int gp = 0; gp < NGP::SQUARE8; ++gp) {
					op.template operator()<8, NGP::SQUARE8>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::TETRA4):
				for (int gp = 0; gp < NGP::TETRA4; ++gp) {
					op.template operator()<4, NGP::TETRA4>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::TETRA10):
				for (int gp = 0; gp < NGP::TETRA10; ++gp) {
					op.template operator()<10, NGP::TETRA10>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID5):
				for (int gp = 0; gp < NGP::PYRAMID5; ++gp) {
					op.template operator()<5, NGP::PYRAMID5>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID13):
				for (int gp = 0; gp < NGP::PYRAMID13; ++gp) {
					op.template operator()<13, NGP::PYRAMID13>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA6):
				for (int gp = 0; gp < NGP::PRISMA6; ++gp) {
					op.template operator()<6, NGP::PRISMA6>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA15):
				for (int gp = 0; gp < NGP::PRISMA15; ++gp) {
					op.template operator()<15, NGP::PRISMA15>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::HEXA8):
				for (int gp = 0; gp < NGP::HEXA8; ++gp) {
					op.template operator()<8, NGP::HEXA8>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::HEXA20):
				for (int gp = 0; gp < NGP::HEXA20; ++gp) {
					op.template operator()<20, NGP::HEXA20>(gp);
				}
			break;
			default: break;
			}
		} else {
			if (Operator::print) printf("\tOP::ELEMENT::%d::GP::%s\n", op.interval, op.name());
			switch (info::mesh->elements->eintervals[op.interval].code) {
			case static_cast<int>(Element::CODE::LINE2):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::LINE2; ++gp) {
						op.template operator()<2, NGP::LINE2>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::LINE3):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::LINE3; ++gp) {
						op.template operator()<3, NGP::LINE3>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE3):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TRIANGLE3; ++gp) {
						op.template operator()<3, NGP::TRIANGLE3>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE6):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TRIANGLE6; ++gp) {
						op.template operator()<6, NGP::TRIANGLE6>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE4):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::SQUARE4; ++gp) {
						op.template operator()<4, NGP::SQUARE4>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE8):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::SQUARE8; ++gp) {
						op.template operator()<8, NGP::SQUARE8>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TETRA4):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TETRA4; ++gp) {
						op.template operator()<4, NGP::TETRA4>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TETRA10):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TETRA10; ++gp) {
						op.template operator()<10, NGP::TETRA10>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID5):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PYRAMID5; ++gp) {
						op.template operator()<5, NGP::PYRAMID5>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID13):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PYRAMID13; ++gp) {
						op.template operator()<13, NGP::PYRAMID13>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA6):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PRISMA6; ++gp) {
						op.template operator()<6, NGP::PRISMA6>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA15):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PRISMA15; ++gp) {
						op.template operator()<15, NGP::PRISMA15>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::HEXA8):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::HEXA8; ++gp) {
						op.template operator()<8, NGP::HEXA8>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::HEXA20):
				for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::HEXA20; ++gp) {
						op.template operator()<20, NGP::HEXA20>(gp);
					}
				}
			break;
			default: break;
			}
		}
	} else {
		if (Operator::print > 2) printf("\tOP::ELEMENT::%d::GP::SKIPPED::%s\n", op.interval, op.name());
	}
}

template <class Operator>
static inline void iterate_elements(Operator op)
{
	if (op.update) {
		if (op.isconst) {
			if (Operator::print) printf("\tOP::ELEMENT::%d::CONSTANT::%s\n", op.interval, op.name());
			op();
		} else {
			if (Operator::print) printf("\tOP::ELEMENT::%d::%s\n", op.interval, op.name());
			for (esint i = info::mesh->elements->eintervals[op.interval].begin; i < info::mesh->elements->eintervals[op.interval].end; ++i, ++op) {
				op();
			}
		}
	} else {
		if (Operator::print > 2) printf("\tOP::ELEMENT::%d::SKIPPED::%s\n", op.interval, op.name());
	}
}

template <class NGP, class Operator>
static inline void iterate_boundary_gps(Operator op, int region)
{
	if (op.update) {
		if (op.isconst) {
			if (Operator::print) printf("\tOP::BOUNDARY::%d::%d::GP::CONSTANT::%s\n", region, op.interval, op.name());
			switch (info::mesh->boundaryRegions[region]->eintervals[op.interval].code) {
			case static_cast<int>(Element::CODE::LINE2):
				for (int gp = 0; gp < NGP::LINE2; ++gp) {
					op.template operator()<2, NGP::LINE2>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::LINE3):
			for (int gp = 0; gp < NGP::LINE3; ++gp) {
				op.template operator()<3, NGP::LINE3>(gp);
			}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE3):
				for (int gp = 0; gp < NGP::TRIANGLE3; ++gp) {
					op.template operator()<3, NGP::TRIANGLE3>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE6):
				for (int gp = 0; gp < NGP::TRIANGLE6; ++gp) {
					op.template operator()<6, NGP::TRIANGLE6>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE4):
				for (int gp = 0; gp < NGP::SQUARE4; ++gp) {
					op.template operator()<4, NGP::SQUARE4>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE8):
				for (int gp = 0; gp < NGP::SQUARE8; ++gp) {
					op.template operator()<8, NGP::SQUARE8>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::TETRA4):
				for (int gp = 0; gp < NGP::TETRA4; ++gp) {
					op.template operator()<4, NGP::TETRA4>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::TETRA10):
				for (int gp = 0; gp < NGP::TETRA10; ++gp) {
					op.template operator()<10, NGP::TETRA10>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID5):
				for (int gp = 0; gp < NGP::PYRAMID5; ++gp) {
					op.template operator()<5, NGP::PYRAMID5>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID13):
				for (int gp = 0; gp < NGP::PYRAMID13; ++gp) {
					op.template operator()<13, NGP::PYRAMID13>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA6):
				for (int gp = 0; gp < NGP::PRISMA6; ++gp) {
					op.template operator()<6, NGP::PRISMA6>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA15):
				for (int gp = 0; gp < NGP::PRISMA15; ++gp) {
					op.template operator()<15, NGP::PRISMA15>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::HEXA8):
				for (int gp = 0; gp < NGP::HEXA8; ++gp) {
					op.template operator()<8, NGP::HEXA8>(gp);
				}
			break;
			case static_cast<int>(Element::CODE::HEXA20):
				for (int gp = 0; gp < NGP::HEXA20; ++gp) {
					op.template operator()<20, NGP::HEXA20>(gp);
				}
			break;
			default: break;
			}
		} else {
			if (Operator::print) printf("\tOP::BOUNDARY::%d::%d::GP::%s\n", region, op.interval, op.name());
			switch (info::mesh->boundaryRegions[region]->eintervals[op.interval].code) {
			case static_cast<int>(Element::CODE::LINE2):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::LINE2; ++gp) {
						op.template operator()<2, NGP::LINE2>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::LINE3):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::LINE3; ++gp) {
						op.template operator()<3, NGP::LINE3>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE3):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TRIANGLE3; ++gp) {
						op.template operator()<3, NGP::TRIANGLE3>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TRIANGLE6):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TRIANGLE6; ++gp) {
						op.template operator()<6, NGP::TRIANGLE6>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE4):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::SQUARE4; ++gp) {
						op.template operator()<4, NGP::SQUARE4>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::SQUARE8):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::SQUARE8; ++gp) {
						op.template operator()<8, NGP::SQUARE8>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TETRA4):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TETRA4; ++gp) {
						op.template operator()<4, NGP::TETRA4>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::TETRA10):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::TETRA10; ++gp) {
						op.template operator()<10, NGP::TETRA10>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID5):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PYRAMID5; ++gp) {
						op.template operator()<5, NGP::PYRAMID5>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PYRAMID13):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PYRAMID13; ++gp) {
						op.template operator()<13, NGP::PYRAMID13>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA6):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PRISMA6; ++gp) {
						op.template operator()<6, NGP::PRISMA6>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::PRISMA15):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::PRISMA15; ++gp) {
						op.template operator()<15, NGP::PRISMA15>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::HEXA8):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::HEXA8; ++gp) {
						op.template operator()<8, NGP::HEXA8>(gp);
					}
				}
			break;
			case static_cast<int>(Element::CODE::HEXA20):
				for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
					for (int gp = 0; gp < NGP::HEXA20; ++gp) {
						op.template operator()<20, NGP::HEXA20>(gp);
					}
				}
			break;
			default: break;
			}
		}
	} else {
		if (Operator::print > 2) printf("\tOP::BOUNDARY::%d::%d::GP::SKIPPED::%s\n", region, op.interval, op.name());
	}
}

template <class Operator>
static inline void iterate_boundary(Operator op, int region)
{
	if (op.update) {
		if (op.isconst) {
			if (Operator::print) printf("\tOP::BOUNDARY::%d::%d::CONSTANT::%s\n", region, op.interval, op.name());
			op();
		} else {
			if (Operator::print) printf("\tOP::BOUNDARY::%d::%d::%s\n", region, op.interval, op.name());
			for (esint i = info::mesh->boundaryRegions[region]->eintervals[op.interval].begin; i < info::mesh->boundaryRegions[region]->eintervals[op.interval].end; ++i, ++op) {
				op();
			}
		}
	} else {
		if (Operator::print > 2) printf("\tOP::BOUNDARY::%d::%d::SKIPPED::%s\n", region, op.interval, op.name());
	}
}


}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_HPP_ */
