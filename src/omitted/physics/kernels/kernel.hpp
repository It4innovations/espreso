
#include "kernel.h"

#include "basis/containers/serializededata.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template<class Instance>
void KernelExecutor<Instance>::nextSubstep()
{
	static_cast<Instance*>(this)->iterator.nextSubstep();
	for (size_t i = 0; i < static_cast<Instance*>(this)->boundaries.size(); ++i) {
		static_cast<Instance*>(this)->boundaries[i].nextSubstep();
	}
	solutionChanged();
}

template<class Instance>
void KernelExecutor<Instance>::solutionChanged()
{
	static_cast<Instance*>(this)->iterator.solutionChanged();
	for (size_t i = 0; i < static_cast<Instance*>(this)->boundaries.size(); ++i) {
		static_cast<Instance*>(this)->boundaries[i].solutionChanged();
	}
}

template<class Instance>
bool KernelExecutor<Instance>::boundaryWithSettings(size_t rindex)
{
	return static_cast<Instance*>(this)->boundaries[rindex].hasSettings();
}

template<class Instance>
void KernelExecutor<Instance>::processElements(const Builder &builder, InstanceFiller &filler)
{
	const ElementsInterval &interval = info::mesh->elements->eintervals[filler.interval];
	auto iterator = static_cast<Instance*>(this)->iterator;
	iterator.initKernelIterator(interval.begin);

	for (esint e = interval.begin; e < interval.end; ++e) {
		iterator.offset = e;
		iterator.element = info::mesh->elements->epointers->datatarray()[e];
		iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

		static_cast<Instance*>(this)->processElement(builder, iterator, filler);
		filler.insert();

		iterator.next({ iterator.element->nodes, 1 });
	}
}

template<class Instance>
void KernelExecutor<Instance>::processBoundary(const Builder &builder, size_t rindex, InstanceFiller &filler)
{
	if (info::mesh->boundaryRegions[rindex]->dimension == 2) {
		processFaces(static_cast<Instance*>(this), builder, rindex, filler);
	}

	if (info::mesh->boundaryRegions[rindex]->dimension == 1) {
		processEdges(static_cast<Instance*>(this), builder, rindex, filler);
	}

	if (info::mesh->boundaryRegions[rindex]->dimension == 0) {
		processNodes(static_cast<Instance*>(this), builder, rindex, filler);
	}
}

template<class Instance>
void KernelExecutor<Instance>::processSolution()
{
	nodesSolution(static_cast<Instance*>(this));
	elementsSolution(static_cast<Instance*>(this));
}

template<class Instance>
template<typename T>
typename std::enable_if<std::is_member_function_pointer<decltype(&T::processNode)>::value>::type
KernelExecutor<Instance>::processNodes(T * t, const Builder &builder, size_t rindex, InstanceFiller &filler)
{
	filler.Ke.resize(0, 0);
	auto iterator = static_cast<Instance*>(this)->boundaries[rindex];
	iterator.initKernelIterator(filler.begin);

	for (esint e = filler.begin; e < filler.end; ++e) {
		static_cast<Instance*>(this)->processNode(builder, iterator, filler);
		filler.insert();

		iterator.next({ 1, 0 });
	}
}

template<class Instance>
template<typename T>
typename std::enable_if<std::is_member_function_pointer<decltype(&T::processEdge)>::value>::type
KernelExecutor<Instance>::processEdges(T * t, const Builder &builder, size_t rindex, InstanceFiller &filler)
{
	auto iterator = static_cast<Instance*>(this)->boundaries[rindex];
	iterator.initKernelIterator(filler.begin);

	for (esint e = filler.begin; e < filler.end; ++e) {
		iterator.element = info::mesh->boundaryRegions[rindex]->epointers->datatarray()[e];

		static_cast<Instance*>(this)->processEdge(builder, iterator, filler);
		filler.insert();

		iterator.next({ iterator.element->nodes, 1 });
	}
}

template<class Instance>
template<typename T>
typename std::enable_if<std::is_member_function_pointer<decltype(&T::processFace)>::value>::type
KernelExecutor<Instance>::processFaces(T * t, const Builder &builder, size_t rindex, InstanceFiller &filler)
{
	auto iterator = static_cast<Instance*>(this)->boundaries[rindex];
	iterator.initKernelIterator(filler.begin);

	for (esint e = filler.begin; e < filler.end; ++e) {
		iterator.element = info::mesh->boundaryRegions[rindex]->epointers->datatarray()[e];

		static_cast<Instance*>(this)->processFace(builder, iterator, filler);
		filler.insert();

		iterator.next({ iterator.element->nodes, 1 });
	}
}

template<class Instance>
template<typename T>
typename std::enable_if<std::is_member_function_pointer<decltype(&T::nodeSolution)>::value>::type
KernelExecutor<Instance>::nodesSolution(T * t)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		auto iterator = static_cast<Instance*>(this)->iterator;
		iterator.initOutputIterator(info::mesh->nodes->distribution[t]);

		for (size_t n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n) {
			static_cast<Instance*>(this)->nodeSolution(iterator);
			iterator.next({ 1, 0 });
		}
	}
	solutionChanged();
}

template<class Instance>
template<typename T>
typename std::enable_if<std::is_member_function_pointer<decltype(&T::elementSolution)>::value>::type
KernelExecutor<Instance>::elementsSolution(T * t)
{
	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		auto iterator = static_cast<Instance*>(this)->iterator;
		iterator.initKernelIterator(info::mesh->elements->distribution.threads[t]);

		for (size_t e = info::mesh->elements->distribution.threads[t]; e < info::mesh->elements->distribution.threads[t + 1]; ++e) {
			iterator.element = info::mesh->elements->epointers->datatarray()[e];
			iterator.material = info::mesh->materials[info::mesh->elements->material->datatarray()[e]];

			static_cast<Instance*>(this)->elementSolution(iterator);
			iterator.next({ iterator.element->nodes, 1 });
		}
	}
}

}

