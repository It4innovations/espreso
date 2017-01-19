
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
void CubeUtils<TElement>::computeInterval(const CubeSettings &settings, size_t cluster[], const Interval &interval, size_t start[], size_t end[])
{
	for (size_t i = 0; i < 3; i++) {
		size_t elements = settings.clusters[i] * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i];
		double istart = interval.start[i] < 0 ? 0 : interval.start[i];
		double iend = interval.end[i] > settings.problemLength[i] ? settings.problemLength[i] : interval.end[i];
		if (istart == iend) {
			if (istart == 0) {
				start[i] = 0;
				end[i] = 1;
			} else {
				end[i] = elements;
				start[i] = end[i] - 1;
			}
		} else {
			double s = (istart / settings.problemLength[i]) * elements;
			double e = (iend / settings.problemLength[i]) * elements;
			start[i] = std::floor(s);
			end[i] = std::ceil(e);
		}

		size_t cStart = cluster[i] * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i];
		size_t cEnd = (cluster[i] + 1) * settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i];

		if (end[i] < cStart || cEnd < start[i]) {
			// will be dropped lated
			return;
		}
		start[i] = start[i] < cStart ? 0 : start[i] - cStart;
		end[i] = end[i] > cEnd ? settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] : end[i] - cStart;

		if (start[i] == end[i]) {
			if (start[i] == 0) {
				start[i] = 0;
				end[i] = 1;
			} else {
				end[i] = elements;
				start[i] = end[i] - 1;
			}
		}
	}
}

template <class TElement>
CubeEdges CubeUtils<TElement>::cubeEdge(const CubeSettings &settings, size_t cluster[], const Interval &interval)
{
	size_t fixed_x = (interval.start[0] == interval.end[0]) ? 1 : 0;
	size_t fixed_y = (interval.start[1] == interval.end[1]) ? 1 : 0;
	size_t fixed_z = (interval.start[2] == interval.end[2]) ? 1 : 0;

	if (fixed_x + fixed_y + fixed_z < 2) {
		return CubeEdges::NONE;
	}

	if (fixed_z) {
		if (interval.start[2]) {
			// Z == 1
			if (fixed_y) {
				if (interval.start[1]) {
					// Z == 1 && Y == 1
					if (cluster[1] == settings.clusters[1] - 1 && cluster[2] == settings.clusters[2] - 1) {
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
				if (interval.start[0]) {
					// Z == 1 && X == 1
					if (cluster[0] == settings.clusters[0] - 1 && cluster[2] == settings.clusters[2] - 1) {
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
				if (interval.start[1]) {
					// Z == 0 && Y == 1
					if (cluster[1] == settings.clusters[1] - 1 && cluster[2] == 0) {
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
				if (interval.start[0]) {
					// Z == 0 && X == 1
					if (cluster[0] == settings.clusters[0] - 1 && cluster[2] == 0) {
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
		if (interval.start[0]) {
			// X == 1
			if (interval.start[1]) {
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
			if (interval.start[1]) {
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

template <class TElement>
CubeFaces CubeUtils<TElement>::cubeFace(const CubeSettings &settings, size_t cluster[], const Interval &interval)
{
	size_t fixed_x = (interval.start[0] == interval.end[0]) ? 1 : 0;
	size_t fixed_y = (interval.start[1] == interval.end[1]) ? 1 : 0;
	size_t fixed_z = (interval.start[2] == interval.end[2]) ? 1 : 0;

	if (fixed_x + fixed_y + fixed_z != 1) {
		return CubeFaces::NONE;
	}

	if (fixed_x) {
		if (interval.start[0] == 0) {
			if (cluster[0] == 0) {
				return CubeFaces::X_0;
			}
		} else {
			if (cluster[0] == settings.clusters[0] - 1) {
				return CubeFaces::X_1;
			}
		}
	}

	if (fixed_y) {
		if (interval.start[1] == 0) {
			if (cluster[1] == 0) {
				return CubeFaces::Y_0;
			}
		} else {
			if (cluster[1] == settings.clusters[1] - 1) {
				return CubeFaces::Y_1;
			}
		}
	}

	if (fixed_z) {
		if (interval.start[2] == 0) {
			if (cluster[2] == 0) {
				return CubeFaces::Z_0;
			}
		} else {
			if (cluster[2] == settings.clusters[2] - 1) {
				return CubeFaces::Z_1;
			}
		}
	}



	return CubeFaces::NONE;
}

}
}
