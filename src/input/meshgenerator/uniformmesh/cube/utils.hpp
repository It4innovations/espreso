
#include "utils.h"

namespace espreso {
namespace input {

template <class TElement>
void CubeUtils<TElement>::globalNodesCount(const CubeSettings &settings, esglobal nodes[])
{
	eslocal cluster[3];
	UniformUtils<TElement>::clusterNodesCount(settings, cluster);
	for (eslocal i = 0; i < 3; i++) {
		nodes[i] = settings.clusters[i] * (cluster[i] - 1) + 1;
	}
}

template <class TElement>
void CubeUtils<TElement>::computeInterval(const CubeSettings &settings, const Interval &interval, size_t start[], size_t end[])
{
	start[0] = std::round((interval.getStart(0) / settings.problemLength[0]) * settings.subdomainsInCluster[0] * settings.elementsInSubdomain[0]);
	start[1] = std::round((interval.getStart(1) / settings.problemLength[1]) * settings.subdomainsInCluster[1] * settings.elementsInSubdomain[1]);
	start[2] = std::round((interval.getStart(2) / settings.problemLength[2]) * settings.subdomainsInCluster[2] * settings.elementsInSubdomain[2]);

	end[0] = std::round((interval.getEnd(0) / settings.problemLength[0]) * settings.subdomainsInCluster[0] * settings.elementsInSubdomain[0]);
	end[1] = std::round((interval.getEnd(1) / settings.problemLength[1]) * settings.subdomainsInCluster[1] * settings.elementsInSubdomain[1]);
	end[2] = std::round((interval.getEnd(2) / settings.problemLength[2]) * settings.subdomainsInCluster[2] * settings.elementsInSubdomain[2]);
}

template <class TElement>
CubeEdges CubeUtils<TElement>::cubeEdge(const CubeSettings &settings, size_t cluster[], size_t start[], size_t end[])
{
	size_t fixed_x = (start[0] == end[0]) ? 1 : 0;
	size_t fixed_y = (start[1] == end[1]) ? 1 : 0;
	size_t fixed_z = (start[2] == end[2]) ? 1 : 0;

	if (fixed_x + fixed_y + fixed_z < 2) {
		ESINFO(GLOBAL_ERROR) << "Incorrect edge interval";
	}

	if (fixed_z) {
		if (start[2]) {
			// Z == 1
			if (fixed_y) {
				if (start[1]) {
					// Z == 1 && Y == 1
					if (cluster[1] == settings.clusters[0] - 1 && cluster[2] == settings.clusters[2] - 1) {
						return CubeEdges::Y_1_Z_1;
					}
				} else {
					// Z == 1 && Y == 0
					if (cluster[1] == 0 && cluster[2] == settings.clusters[2] - 1) {
						return CubeEdges::Y_0_Z_1;
					}
				}
			}
			if (fixed_x) {
				if (start[0]) {
					// Z == 1 && X == 1
					if (cluster[0] == settings.clusters[1] - 1 && cluster[2] == settings.clusters[2] - 1) {
						return CubeEdges::X_1_Z_1;
					}
				} else {
					// Z == 1 && X == 0
					if (cluster[0] == 0 && cluster[2] == settings.clusters[2] - 1) {
						return CubeEdges::X_0_Z_1;
					}
				}
			}
		} else {
			// Z == 0
			if (fixed_y) {
				if (start[1]) {
					// Z == 0 && Y == 1
					if (cluster[1] == settings.clusters[0] - 1 && cluster[2] == 0) {
						return CubeEdges::Y_1_Z_0;
					}
				} else {
					// Z == 0 && Y == 0
					if (cluster[1] == 0 && cluster[2] == 0) {
						return CubeEdges::Y_0_Z_0;
					}
				}
			}
			if (fixed_x) {
				if (start[0]) {
					// Z == 0 && X == 1
					if (cluster[0] == settings.clusters[1] - 1 && cluster[2] == 0) {
						return CubeEdges::X_1_Z_0;
					}
				} else {
					// Z == 0 && X == 0
					if (cluster[0] == 0 && cluster[2] == 0) {
						return CubeEdges::X_0_Z_0;
					}
				}
			}

		}
	} else {
		if (start[0]) {
			// X == 1
			if (start[1]) {
				// X == 1 && Y == 1
				if (cluster[0] == settings.clusters[0] - 1 && cluster[1] == settings.clusters[1] - 1) {
					return CubeEdges::X_1_Y_1;
				}
			} else {
				// X == 1 && Y == 0
				if (cluster[0] == settings.clusters[0] - 1 && cluster[1] == 0) {
					return CubeEdges::X_1_Y_0;
				}
			}
		} else {
			// X == 0
			if (start[1]) {
				// X == 0 && Y == 1
				if (cluster[0] == 0 && cluster[1] == settings.clusters[1] - 1) {
					return CubeEdges::X_0_Y_1;
				}
			} else {
				// X == 0 && Y == 0
				if (cluster[0] == 0 && cluster[1] == 0) {
					return CubeEdges::X_0_Y_0;
				}
			}
		}
	}

	return CubeEdges::NONE;
}

}
}
