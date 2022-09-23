
#include "builder.utils.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "wrappers/mpi/communication.h"

namespace espreso {
namespace builder {

void clip(LinkedNodes &nodes, MergedElements &elements)
{
//	if (!info::ecf->input.clipping_box.apply) {
//		return;
//	}
//
//	eslog::info(" == CLIPPED BOX           <%9f %9f> <%9f %9f> <%9f %9f> == \n",
//			info::ecf->input.clipping_box.min[0], info::ecf->input.clipping_box.max[0],
//			info::ecf->input.clipping_box.min[1], info::ecf->input.clipping_box.max[1],
//			info::ecf->input.clipping_box.min[2], info::ecf->input.clipping_box.max[2]);
//	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
//
//	esint esize = 0, last = 0;
//	for (esint e = 0, eoffset = 0; e < elements.size; eoffset += Element::encode(elements.etype[e++]).nodes) {
//		bool isIn = true;
//		PolyElement poly(elements.etype[e], elements.enodes.data() + eoffset);
//		for (esint n = 0; n < Element::encode(elements.etype[e]).nodes; ++n) {
//			if (poly.isNode(n)) {
//				const Point &p = nodes.coordinates[elements.enodes[n + eoffset] - nodes.offset];
//				if (
//						p.x < info::ecf->input.clipping_box.min[0] || info::ecf->input.clipping_box.max[0] < p.x ||
//						p.y < info::ecf->input.clipping_box.min[1] || info::ecf->input.clipping_box.max[1] < p.y ||
//						p.z < info::ecf->input.clipping_box.min[2] || info::ecf->input.clipping_box.max[2] < p.z) {
//					isIn = false;
//				}
//			}
//		}
//		if (isIn) {
//			for (esint n = 0; n < Element::encode(elements.etype[e]).nodes; ++n) {
//				elements.enodes[last++] = elements.enodes[n + eoffset];
//			}
//			elements.etype[esize++] = elements.etype[e];
//		}
//	}
//	elements.etype.resize(esize);
//	elements.enodes.resize(last);
//	elements.size = esize;
//	elements.offset = esize;
//	elements.total = Communication::exscan(elements.offset);
}

}
}
