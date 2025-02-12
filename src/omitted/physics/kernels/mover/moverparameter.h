
#ifndef SRC_PHYSICS_KERNELS_MOVER_MOVERPARAMETER_H_
#define SRC_PHYSICS_KERNELS_MOVER_MOVERPARAMETER_H_

#include "mover.h"
#include "basis/containers/serializededata.h"
#include "esinfo/envinfo.h"
#include "mesh/store/nameddata.h"

namespace espreso {

class MoverParameter {
public:
    enum Properties {
        UNKNOWN,
        ALLOW_CONSTANT,
    };

    int dimension;
    double defaultValue;
    const char* name;

    double *data;
    MoverIncrement increment;

    MoverParameter(int dimension, double defaultValue, const char* name = NULL)
    : dimension(dimension), defaultValue(defaultValue), name(name), data(NULL) {}
    virtual ~MoverParameter() {}

    virtual void clear() =0;
    virtual void initKernelIterator(esint offset) =0;
    virtual void initOutputIterator(esint offset) =0;
    virtual void updatingFromInput(MoverInstance *iterator) =0;
    virtual void updatingFromOutput(MoverInstance *iterator) =0;

protected:
    template <class TParameter>
    void setDataPointer(TParameter &parameter, esint offset);

    template <class TInput, class TKernel, class TECFExpr>
    void setInput(TInput &input, TKernel &kernel, TECFExpr &ecf, Evaluator::Params &params, serializededata<esint, esint> *nodes, Properties properties)
    {
        input.ecf = &ecf;
        input.params = &params;

        kernel.dimension = dimension;
        kernel.nodes = nodes;
        kernel.isConst = (properties & Properties::ALLOW_CONSTANT) && input.isConst();
        if (kernel.isConst) {
            kernel.values = new serializededata<esint, double>(dimension, tarray<double>(info::env::OMP_NUM_THREADS, dimension * info::env::OMP_NUM_THREADS * 64, defaultValue));
        } else {
            kernel.values = new serializededata<esint, double>(dimension, nodes->datatarray().distribution(), defaultValue);
        }
        data = kernel.values->datatarray().data();
    }

    template <class TOutput>
    void setOutput(TOutput &output, NamedData::DataType datatype, const char* name, bool guard, step::TYPE restriction);
};

template <class TInput, class TOutput>
class MoverFullParameter: public MoverParameter {
public:
    TInput input;
    ElementNodeValues kernel;
    TOutput output;

    MoverFullParameter(int dimension, double defaultValue, const char* name = NULL)
    : MoverParameter(dimension, defaultValue, name) {}

    void clear()
    {
        if (kernel.values) { delete kernel.values; }
    }

    void initKernelIterator(esint offset)
    {
        MoverParameter::setDataPointer(kernel, offset);
    }

    void initOutputIterator(esint offset)
    {
        MoverParameter::setDataPointer(output, offset);
    }

    void initFrom(decltype(*input.ecf) &ecf)
    {
        TInput initialization;
        initialization.ecf = &ecf;
        initialization.params = input.params;

        if (kernel.values) {
            Move<TInput, ElementNodeValues>(initialization, kernel)();
        }
        if (output.data) {
            Move<ElementNodeValues, TOutput>(kernel, output)();
        }
    }

    void updatingFromInput(MoverInstance *iterator)
    {
        if (kernel.values) {
            iterator->afterSolutionChange(input, kernel);
        }
        if (output.data) {
            iterator->afterSolutionChange(kernel, output);
        }
    }

    void updatingFromOutput(MoverInstance *iterator)
    {
        if (output.data && kernel.values) {
            iterator->afterSolutionChange(output, kernel);
        }
    }

    void setInput(
            decltype(*input.ecf) &ecf,
            Evaluator::Params &params,
            serializededata<esint, esint> *nodes,
            Properties properties = Properties::UNKNOWN)
    {
        MoverParameter::setInput(input, kernel, ecf, params, nodes, properties);
    }

    void setOutput(
            NamedData::DataType datatype,
            const char* name,
            bool guard = true,
            step::TYPE restriction = step::TYPE::TIME)
    {
        MoverParameter::setOutput<TOutput>(output, datatype, name, guard, restriction);
    }
};

class MoverCoordinates: public MoverParameter {
public:
    InputCoordinates input;
    ElementNodeValues kernel;

    MoverCoordinates(int dimension, const char* name = NULL)
    : MoverParameter(dimension, 0, name) {}

    void clear();
    void initKernelIterator(esint offset);
    void initOutputIterator(esint offset);
    void updatingFromInput(MoverInstance *iterator);
    void updatingFromOutput(MoverInstance *iterator);
    void set(serializededata<esint, Point>* points, serializededata<esint, esint> *nodes);
};

template <class TInput>
class MoverInputParameter: public MoverParameter {
public:
    TInput input;
    ElementNodeValues kernel;

    MoverInputParameter(int dimension, double defaultValue, const char* name = NULL)
    : MoverParameter(dimension, defaultValue, name) {}

    void clear()
    {
        if (kernel.values) { delete kernel.values; }
    }

    void initKernelIterator(esint offset)
    {
        MoverParameter::setDataPointer(kernel, offset);
    }

    void initOutputIterator(esint offset)
    {
        data = NULL;
    }

    void updatingFromInput(MoverInstance *iterator)
    {
        if (kernel.values) {
            iterator->afterSolutionChange(input, kernel);
        }
    }

    void updatingFromOutput(MoverInstance *iterator)
    {

    }

    template <class TECFInput>
    void setInput(
            TECFInput &ecf,
            Evaluator::Params &params,
            serializededata<esint, esint> *nodes,
            Properties properties = Properties::UNKNOWN)
    {
        MoverParameter::setInput(input, kernel, ecf, params, nodes, properties);
    }

    template <class TECFInput>
    void setInput(
            std::map<std::string, TECFInput> &ecf,
            const std::string &name,
            Evaluator::Params &params,
            serializededata<esint, esint> *nodes,
            Properties properties = Properties::UNKNOWN)
    {
        if (ecf.find(name) != ecf.end()) {
            MoverParameter::setInput(input, kernel, ecf.find(name)->second, params, nodes, properties);
        }
    }
};

template <class ...TArgs>
class MoverInputStructParameter {
public:

    void clear()
    {
        for (size_t i = 0; i < _params.size(); ++i) {
            _params[i]->clear();
        }
    }

    void initKernelIterator(esint offset)
    {
        for (size_t i = 0; i < _params.size(); ++i) {
            _params[i]->initKernelIterator(offset);
        }
    }

    void initOutputIterator(esint offset)
    {
        for (size_t i = 0; i < _params.size(); ++i) {
            _params[i]->initOutputIterator(offset);
        }
    }

    void updatingFromInput(MoverInstance *iterator)
    {
        for (size_t i = 0; i < _params.size(); ++i) {
            _params[i]->updatingFromInput(iterator);
        }
    }

    void updatingFromOutput(MoverInstance *iterator)
    {
        for (size_t i = 0; i < _params.size(); ++i) {
            _params[i]->updatingFromOutput(iterator);
        }
    }

protected:
    std::vector<MoverParameter*> _params;
};

class MoverInputHarmonicParameter: public MoverInputStructParameter<InputHarmonicMagnitudeExpressionVector, InputHarmonicPhaseExpressionVector> {
public:
    MoverInputParameter<InputHarmonicMagnitudeExpressionVector> magnitude;
    MoverInputParameter<InputHarmonicPhaseExpressionVector> phase;

    MoverInputHarmonicParameter(int dimension, double defaultValue, const char* name = NULL)
    : magnitude(dimension, defaultValue, name), phase(dimension, defaultValue, name)
    {
        _params.push_back(&magnitude);
        _params.push_back(&phase);
    }

    void setInput(
            std::map<std::string, ECFHarmonicExpressionVector> &ecf,
            const std::string &name,
            Evaluator::Params &params,
            serializededata<esint, esint> *nodes,
            MoverParameter::Properties properties = MoverParameter::Properties::UNKNOWN)
    {
        if (ecf.find(name) != ecf.end()) {
            magnitude.setInput(ecf.find(name)->second, params, nodes, properties);
            phase.setInput(ecf.find(name)->second, params, nodes, properties);
        }
    }
};

template <class TOutput>
class MoverOutputParameter: public MoverParameter {
public:
    TOutput output;

    MoverOutputParameter(int dimension, double defaultValue, const char* name = NULL)
    : MoverParameter(dimension, defaultValue, name) {}

    void clear()
    {

    }

    void initKernelIterator(esint offset);
    void initOutputIterator(esint offset);

    void updatingFromInput(MoverInstance *iterator)
    {

    }

    void updatingFromOutput(MoverInstance *iterator)
    {

    }

    void setOutput(
            NamedData::DataType datatype,
            const char* name,
            bool guard = true,
            step::TYPE restriction = step::TYPE::TIME | step::TYPE::FREQUENCY | step::TYPE::FTT)
    {
        MoverParameter::setOutput<TOutput>(output, datatype, name, guard, restriction);
    }
};

template<class TInput>
class MoverReducedInputParameter: public MoverParameter {
public:
    TInput input;
    ElementNodeValues kernel;

    MoverReducedInputParameter(int dimension, double defaultValue, const char* name = NULL)
    : MoverParameter(dimension, defaultValue, name) {}

    void clear()
    {
        if (kernel.values) { delete kernel.values; }
    }

    void initKernelIterator(esint offset)
    {
        MoverParameter::setDataPointer(kernel, offset);
    }

    void initOutputIterator(esint eoffset)
    {
        data = NULL;
    }

    void updatingFromInput(MoverInstance *iterator)
    {
        iterator->afterSolutionChange(input, kernel);
    }

    void updatingFromOutput(MoverInstance *iterator)
    {

    }

    void setInput(TInput input, serializededata<esint, esint> *nodes, Properties properties = Properties::UNKNOWN)
    {
        this->input = input;

        kernel.dimension = dimension;
        kernel.nodes = nodes;
        kernel.values = new serializededata<esint, double>(dimension, nodes->datatarray().distribution());
    }
};

template<class TInput, class TOutput>
class MoverBoundaryParameter: public MoverParameter {
public:
    TInput input;
    ElementNodeValues kernel;
    TOutput output;

    MoverBoundaryParameter(int dimension, double defaultValue, const char* name = NULL)
    : MoverParameter(dimension, defaultValue, name) {}

    void clear()
    {
        if (kernel.values) { delete kernel.values; }
    }

    void initKernelIterator(esint offset)
    {
        MoverParameter::setDataPointer(kernel, offset);
    }

    void initOutputIterator(esint offset)
    {
        MoverParameter::setDataPointer(output, offset);
    }

    void updatingFromInput(MoverInstance *iterator)
    {
        iterator->afterSolutionChange(input, kernel);
        iterator->afterSolutionChange(kernel, output);
    }

    void updatingFromOutput(MoverInstance *iterator)
    {

    }

    template <class TECFInput>
    void setInput(
            TECFInput &ecf,
            Evaluator::Params &params,
            serializededata<esint, esint> *nodes,
            Properties properties = Properties::UNKNOWN)
    {
        MoverParameter::setInput(input, kernel, ecf, params, nodes, properties);
    }

    template <class TECFInput>
    void setInput(
            std::map<std::string, TECFInput> &ecf,
            const std::string &name,
            Evaluator::Params &params,
            serializededata<esint, esint> *nodes,
            Properties properties = Properties::UNKNOWN)
    {
        if (ecf.find(name) != ecf.end()) {
            MoverParameter::setInput(input, kernel, ecf.find(name)->second, params, nodes, properties);
        }
    }

    void setOutput(
            NamedData::DataType datatype,
            const char* name,
            bool guard = true,
            step::TYPE restriction = step::TYPE::TIME)
    {
        MoverParameter::setOutput<TOutput>(output, datatype, name, guard, restriction);
    }
};

}



#endif /* SRC_PHYSICS_KERNELS_MOVER_MOVERPARAMETER_H_ */
