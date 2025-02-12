
#include "moverparameter.h"
#include "mover.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"

using namespace espreso;

MoverInstance::MoverInstance()
: element(NULL)
{

}

MoverInstance::MoverInstance(const MoverInstance &other)
: element(NULL)
{
    uintptr_t _this = reinterpret_cast<uintptr_t>(this);
    uintptr_t _other = reinterpret_cast<uintptr_t>(&other);
    inputs = other.inputs;
    for (size_t i = 0; i < other.inputs.size(); i++) {
        uintptr_t _data = reinterpret_cast<uintptr_t>(other.inputs[i]);
        inputs[i] = reinterpret_cast<MoverParameter*>(_this + (_data - _other));
    }
    outputs = other.outputs;
    for (size_t i = 0; i < other.outputs.size(); i++) {
        uintptr_t _data = reinterpret_cast<uintptr_t>(other.outputs[i]);
        outputs[i] = reinterpret_cast<MoverParameter*>(_this + (_data - _other));
    }
}

MoverInstance::~MoverInstance()
{
    for (size_t i = 0; i < stepMovers.size(); i++) {
        delete stepMovers[i];
    }
    for (size_t i = 0; i < solutionMovers.size(); i++) {
        delete solutionMovers[i];
    }
}

void MoverInstance::registerSolution(MoverParameter &parameter, const MoverIncrement &increment)
{
    parameter.increment = increment;
    inputs.push_back(&parameter);
    parameter.updatingFromOutput(this);
}

void MoverInstance::registerInput(MoverParameter &parameter, const MoverIncrement &increment)
{
    parameter.increment = increment;
    inputs.push_back(&parameter);
    parameter.updatingFromInput(this);
}

void MoverInstance::registerOutput(MoverParameter &parameter, const MoverIncrement &increment)
{
    parameter.increment = increment;
    outputs.push_back(&parameter);
}

void MoverInstance::clear()
{
    for (size_t i = 0; i < inputs.size(); i++) {
        inputs[i]->clear();
    }
    for (size_t i = 0; i < outputs.size(); i++) {
        outputs[i]->clear();
    }
}

void MoverInstance::initKernelIterator(esint offset)
{
    for (size_t i = 0; i < inputs.size(); i++) {
        inputs[i]->initKernelIterator(offset);
    }
    for (size_t i = 0; i < outputs.size(); i++) {
        outputs[i]->initKernelIterator(offset);
    }
}

void MoverInstance::initOutputIterator(esint offset)
{
    for (size_t i = 0; i < inputs.size(); i++) {
        inputs[i]->initOutputIterator(offset);
    }
    for (size_t i = 0; i < outputs.size(); i++) {
        outputs[i]->initOutputIterator(offset);
    }
}

void MoverInstance::next(const MoverIncrement &counters)
{
    for (size_t i = 0; i < inputs.size(); i++) {
        if (inputs[i]->data) {
            inputs[i]->data += inputs[i]->increment * counters;
        }
    }
    for (size_t i = 0; i < outputs.size(); i++) {
        if (outputs[i]->data) {
            outputs[i]->data += outputs[i]->increment * counters;
        }
    }
}

void MoverInstance::nextSubstep()
{
    nodeparams.time(step::time.current);
    nodeparams.freq(step::frequency.current);
    kernelparams.time(step::time.current);
    kernelparams.freq(step::frequency.current);
    for (size_t i = 0; i < stepMovers.size(); i++) {
        stepMovers[i]->operator()();
    }
}

void MoverInstance::solutionChanged()
{
    nodeparams.time(step::time.current);
    nodeparams.freq(step::frequency.current);
    kernelparams.time(step::time.current);
    kernelparams.freq(step::frequency.current);
    for (size_t i = 0; i < solutionMovers.size(); i++) {
        solutionMovers[i]->operator()();
    }
}

void ElementNodeValues::clear()
{
    if (values) delete values;
}


Move<OutputNodes, ElementNodeValues>::Move(const OutputNodes &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{

}

void Move<OutputNodes, ElementNodeValues>::operator()()
{
    #pragma omp parallel for
    for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
        auto i = to.values->datatarray().begin(t);
        for (auto n = to.nodes->datatarray().cbegin(t); n != to.nodes->datatarray().cend(t); ++n) {
            for (int d = 0; d < from.data->dimension; ++d, ++i) {
                *i = from.data->data[*n * from.data->dimension + d];
            }
        }
    }
}

Move<OutputElements, ElementNodeValues>::Move(const OutputElements &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{

}

void Move<OutputElements, ElementNodeValues>::operator()()
{
//    eslog::internalFailure("call empty function.\n");
}

Move<InputCoordinates, ElementNodeValues>::Move(const InputCoordinates &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{

}

void Move<InputCoordinates, ElementNodeValues>::operator()()
{
    if (from.dimension == 2) {
        #pragma omp parallel for
        for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
            auto c = to.values->begin(t);
            for (auto n = to.nodes->datatarray().cbegin(t); n != to.nodes->datatarray().cend(t); ++n, ++c) {
                c->at(0) = from.values->datatarray()[*n].x;
                c->at(1) = from.values->datatarray()[*n].y;
            }
        }
    }
    if (from.dimension == 3) {
        #pragma omp parallel for
        for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
            auto c = to.values->begin(t);
            for (auto n = to.nodes->datatarray().cbegin(t); n != to.nodes->datatarray().cend(t); ++n, ++c) {
                c->at(0) = from.values->datatarray()[*n].x;
                c->at(1) = from.values->datatarray()[*n].y;
                c->at(2) = from.values->datatarray()[*n].z;
            }
        }
    }
}

Move<InputExpression, ElementNodeValues>::Move(const InputExpression &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputExpression, ElementNodeValues>::operator()()
{
    if (from.ecf == NULL || moved) {
        return;
    }
    from.ecf->evaluator->evalVector(
            to.values->datatarray().size(), *from.params,
            to.values->datatarray().data());
}

Move<InputExpressionVector, ElementNodeValues>::Move(const InputExpressionVector &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputExpressionVector, ElementNodeValues>::operator()()
{
    if (from.ecf == NULL || moved) {
        return;
    }
    from.ecf->x.evaluator->evalVector(
            to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
            to.values->datatarray().data() + 0);
    from.ecf->y.evaluator->evalVector(
            to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
            to.values->datatarray().data() + 1);
    if (to.dimension == 3) {
        from.ecf->z.evaluator->evalVector(
                to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                to.values->datatarray().data() + 2);
    }
}

Move<InputExpressionOptionalVector, ElementNodeValues>::Move(const InputExpressionOptionalVector &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputExpressionOptionalVector, ElementNodeValues>::operator()()
{
    if (from.ecf == NULL || moved) {
        return;
    }
    if (from.ecf->all.value.size()) {
        for (int d = 0; d < to.dimension; ++d) {
            from.ecf->all.evaluator->evalVector(
                    to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                    to.values->datatarray().data() + d);
        }
    } else {
        if (from.ecf->x.value.size()) {
            from.ecf->x.evaluator->evalVector(
                    to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                    to.values->datatarray().data() + 0);
        }
        if (from.ecf->y.value.size()) {
            from.ecf->y.evaluator->evalVector(
                    to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                    to.values->datatarray().data() + 1);
        }
        if (to.dimension == 3 && from.ecf->z.value.size()) {
            from.ecf->z.evaluator->evalVector(
                    to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                    to.values->datatarray().data() + 2);
        }
    }
}

Move<InputHarmonicMagnitudeExpression, ElementNodeValues>::Move(const InputHarmonicMagnitudeExpression &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputHarmonicMagnitudeExpression, ElementNodeValues>::operator()()
{
    if (from.ecf == NULL || moved) {
        return;
    }
    from.ecf->magnitude.evaluator->evalVector(
            to.values->datatarray().size(), *from.params,
            to.values->datatarray().data());
}

Move<InputHarmonicPhaseExpression, ElementNodeValues>::Move(const InputHarmonicPhaseExpression &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputHarmonicPhaseExpression, ElementNodeValues>::operator()()
{
    if (from.ecf == NULL || moved) {
        return;
    }
    from.ecf->phase.evaluator->evalVector(
            to.values->datatarray().size(), *from.params,
            to.values->datatarray().data());
}

Move<InputHarmonicMagnitudeExpressionVector, ElementNodeValues>::Move(const InputHarmonicMagnitudeExpressionVector &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputHarmonicMagnitudeExpressionVector, ElementNodeValues>::operator()()
{
    if (from.ecf == NULL || moved) {
        return;
    }
    from.ecf->magnitude.x.evaluator->evalVector(
            to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
            to.values->datatarray().data() + 0);
    from.ecf->magnitude.y.evaluator->evalVector(
            to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
            to.values->datatarray().data() + 1);
    if (to.dimension == 3) {
        from.ecf->magnitude.z.evaluator->evalVector(
                to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                to.values->datatarray().data() + 2);
    }
}

Move<InputHarmonicPhaseExpressionVector, ElementNodeValues>::Move(const InputHarmonicPhaseExpressionVector &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputHarmonicPhaseExpressionVector, ElementNodeValues>::operator()()
{
    if (from.ecf == NULL || moved) {
        return;
    }
    from.ecf->phase.x.evaluator->evalVector(
            to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
            to.values->datatarray().data() + 0);
    from.ecf->phase.y.evaluator->evalVector(
            to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
            to.values->datatarray().data() + 1);
    if (to.dimension == 3) {
        from.ecf->phase.z.evaluator->evalVector(
                to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                to.values->datatarray().data() + 2);
    }
}

Move<InputExpressionMap, OutputNodes>::Move(const InputExpressionMap &from, const OutputNodes &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputExpressionMap, OutputNodes>::operator()()
{
    if (moved) { return; }

    #pragma omp parallel for
    for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
        for (auto it = from.ecf->begin(); it != from.ecf->end(); ++it) {
            serializededata<esint, esint>* nodes = info::mesh->bregion(it->first)->nodes;
            it->second.evaluator->evalSelectedDense(
                    nodes->datatarray().size(t),
                    nodes->datatarray().begin(t),
                    *from.params,
                    to.data->data.data());
        }
    }
}

Move<InputExpressionMap, ElementNodeValues>::Move(const InputExpressionMap &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputExpressionMap, ElementNodeValues>::operator()()
{
    if (moved) { return; }

    #pragma omp parallel for
    for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
        for (auto it = from.ecf->begin(); it != from.ecf->end(); ++it) {
            ElementsRegionStore *region = info::mesh->eregion(it->first);
            if (to.isConst) {
                it->second.evaluator->evalVector(
                        to.values->datatarray().size(), *from.params,
                        to.values->datatarray().data());
            } else {
                it->second.evaluator->evalFiltered(
                        region->elements->datatarray().size(t),
                        region->elements->datatarray().begin(t),
                        to.nodes->boundarytarray().cbegin(),
                        *from.params,
                        to.values->datatarray().data());
            }
        }
    }
}

Move<InputExpressionVectorMap, ElementNodeValues>::Move(const InputExpressionVectorMap &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputExpressionVectorMap, ElementNodeValues>::operator()()
{
    if (moved) { return; }

    #pragma omp parallel for
    for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
        for (auto it = from.ecf->begin(); it != from.ecf->end(); ++it) {
            ElementsRegionStore *region = info::mesh->eregion(it->first);
            if (to.isConst) {
                it->second.x.evaluator->evalVector(
                        to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                        to.values->datatarray().data() + 0);
                it->second.y.evaluator->evalVector(
                        to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                        to.values->datatarray().data() + 1);
                if (to.dimension == 3) {
                    it->second.z.evaluator->evalVector(
                            to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                            to.values->datatarray().data() + 2);
                }
            } else {
                it->second.x.evaluator->evalFiltered(
                        region->elements->datatarray().size(t), from.params->ncoords(),
                        region->elements->datatarray().begin(t),
                        to.nodes->boundarytarray().cbegin(),
                        *from.params,
                        to.values->datatarray().data() + 0
                );
                it->second.y.evaluator->evalFiltered(
                        region->elements->datatarray().size(t), from.params->ncoords(),
                        region->elements->datatarray().begin(t),
                        to.nodes->boundarytarray().cbegin(),
                        *from.params,
                        to.values->datatarray().data() + 1
                );
                if (from.params->ncoords() == 3) {
                    it->second.z.evaluator->evalFiltered(
                            region->elements->datatarray().size(t), from.params->ncoords(),
                            region->elements->datatarray().begin(t),
                            to.nodes->boundarytarray().cbegin(),
                            *from.params,
                            to.values->datatarray().data() + 2
                    );
                }
            }
        }
    }
}

Move<InputExpressionOptionalVectorMap, ElementNodeValues>::Move(const InputExpressionOptionalVectorMap &from, const ElementNodeValues &to)
: from(from), to(to), moved(false)
{
    if (from.isConst()) { now(); moved = true; }
}

void Move<InputExpressionOptionalVectorMap, ElementNodeValues>::operator()()
{
    if (moved) { return; }

    #pragma omp parallel for
    for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
        for (auto it = from.ecf->begin(); it != from.ecf->end(); ++it) {
            ElementsRegionStore *region = info::mesh->eregion(it->first);
            if (to.isConst) {
                if (it->second.all.value.size()) {
                    for (int d = 0; d < to.dimension; ++d) {
                        it->second.all.evaluator->evalVector(
                                to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                                to.values->datatarray().data() + d);
                    }
                } else {
                    if (it->second.x.value.size()) {
                        it->second.x.evaluator->evalVector(
                                to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                                to.values->datatarray().data() + 0);
                    }
                    if (it->second.y.value.size()) {
                        it->second.y.evaluator->evalVector(
                                to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                                to.values->datatarray().data() + 1);
                    }
                    if (to.dimension == 3 && it->second.z.value.size()) {
                        it->second.z.evaluator->evalVector(
                                to.values->datatarray().size() / to.dimension, to.dimension, *from.params,
                                to.values->datatarray().data() + 2);
                    }
                }
            } else {
                if (it->second.all.value.size()) {
                    for (int d = 0; d < from.params->ncoords(); ++d) {
                        it->second.all.evaluator->evalFiltered(
                                region->elements->datatarray().size(t), from.params->ncoords(),
                                region->elements->datatarray().begin(t),
                                to.nodes->boundarytarray().cbegin(),
                                *from.params,
                                to.values->datatarray().data() + d
                        );
                    }
                } else {
                    if (it->second.x.value.size()) {
                        it->second.x.evaluator->evalFiltered(
                                region->elements->datatarray().size(t), from.params->ncoords(),
                                region->elements->datatarray().begin(t),
                                to.nodes->boundarytarray().cbegin(),
                                *from.params,
                                to.values->datatarray().data() + 0
                        );
                    }
                    if (it->second.y.value.size()) {
                        it->second.y.evaluator->evalFiltered(
                                region->elements->datatarray().size(t), from.params->ncoords(),
                                region->elements->datatarray().begin(t),
                                to.nodes->boundarytarray().cbegin(),
                                *from.params,
                                to.values->datatarray().data() + 1
                        );
                    }
                    if (from.params->ncoords() == 3 && it->second.z.value.size()) {
                        it->second.z.evaluator->evalFiltered(
                                region->elements->datatarray().size(t), from.params->ncoords(),
                                region->elements->datatarray().begin(t),
                                to.nodes->boundarytarray().cbegin(),
                                *from.params,
                                to.values->datatarray().data() + 2
                        );
                    }
                }
            }
        }
    }
}

Move<ElementNodeValues, OutputNodes>::Move(const ElementNodeValues &from, const OutputNodes &to)
: from(from), to(to), moved(false)
{

}

void Move<ElementNodeValues, OutputNodes>::operator()()
{
    if (from.values == NULL) {
        return;
    }
    if (to.data == NULL) {
        return;
    }
    if (from.isConst) {
        auto i = from.values->datatarray().cbegin();
        for (size_t n = 0; n < to.data->data.size() / to.data->dimension; n++) {
            for (int d = 0; d < to.data->dimension; d++) {
                to.data->data[n * to.data->dimension + d] = i[d];
            }
        }
        return;
    }

    std::fill(to.data->data.begin(), to.data->data.end(), 0);
    auto i = from.values->datatarray().cbegin();
    for (auto n = from.nodes->datatarray().cbegin(); n != from.nodes->datatarray().cend(); ++n, ++i) {
        for (int d = 0; d < to.data->dimension; d++) {
            to.data->data[*n * to.data->dimension + d] += *i;
        }
    }

    auto &nelements = info::mesh->nodes->elements->boundarytarray();
    for (size_t i = 0; i < to.data->data.size() / to.data->dimension; i++) {
        for (int d = 0; d < to.data->dimension; d++) {
            to.data->data[i * to.data->dimension + d] /= nelements[i + 1] - nelements[i];
        }
    }

    std::vector<std::vector<double> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());

    auto nranks = info::mesh->nodes->ranks->begin();
    for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
        if (nranks->size() > 1) {
            esint noffset = 0;
            for (auto r = nranks->begin(); r != nranks->end(); ++r) {
                if (*r != info::mpi::rank) {
                    while (info::mesh->neighbors[noffset] < *r) {
                        ++noffset;
                    }
                    for (int d = 0; d < to.data->dimension; d++) {
                        sBuffer[noffset].push_back(to.data->data[n * to.data->dimension + d]);
                    }
                }
            }
        }
    }

    for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
        rBuffer[n].resize(sBuffer[n].size());
    }

    if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
        eslog::internalFailure("exchange diagonal values.\n");
    }

    nranks = info::mesh->nodes->ranks->begin();
    std::vector<esint> nindex(info::mesh->neighbors.size());
    for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
        if (nranks->size() > 1) {
            esint noffset = 0;
            for (auto r = nranks->begin(); r != nranks->end(); ++r) {
                if (*r != info::mpi::rank) {
                    while (info::mesh->neighbors[noffset] < *r) {
                        ++noffset;
                    }
                    for (int d = 0; d < to.data->dimension; d++) {
                        to.data->data[n * to.data->dimension + d] += rBuffer[noffset][nindex[noffset]++];
                    }
                }
            }
        }
    }
}

Move<ElementNodeValues, OutputElements>::Move(const ElementNodeValues &from, const OutputElements &to)
: from(from), to(to), moved(false)
{

}

void Move<ElementNodeValues, OutputElements>::operator()()
{
    if (from.values == NULL || to.data == NULL) {
        return;
    }
    if (from.isConst) {
        auto i = from.values->datatarray().cbegin();
        for (size_t n = 0; n < to.data->data.size() / to.data->dimension; n++) {
            for (int d = 0; d < to.data->dimension; d++) {
                to.data->data[n * to.data->dimension + d] = i[d];
            }
        }
        return;
    }

    auto elements = info::mesh->elements;

    #pragma omp parallel for
    for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
        size_t noffset = elements->nodes->datatarray().distribution()[t];
        size_t eoffset = elements->distribution.threads[t];
        for (auto enodes = elements->nodes->cbegin(t); enodes != elements->nodes->cend(t); ++enodes, ++eoffset) {
            for (int d = 0; d < to.data->dimension; d++) {
                double sum = 0;
                for (auto n = enodes->begin(); n != enodes->end(); ++n, ++noffset) {
                    sum += from.values->datatarray()[to.data->dimension * noffset + d];
                }
                to.data->data[to.data->dimension * eoffset + d] = sum / enodes->size();
                noffset -= enodes->size();
            }
            noffset += enodes->size();
        }
    }
}
