
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_

#include "subkernel.h"
#include "basis/evaluator/evaluator.h"
#include "config/holders/expression.h"

#include "omp.h"
#include <functional>

namespace espreso {

struct ExternalExpression {
    const char* name() const { return expression->value.c_str(); }

    ECFExpression *expression;

    ExternalExpression()
    : expression(nullptr)
    {

    }

    void activate(ECFExpression &expression)
    {
        this->expression = &expression;
    }
};

struct ExternalExpressionVector {
    const char* name() const { return "ExternalExpressionVector"; }

    ECFExpressionVector *expression;

    ExternalExpressionVector()
    : expression(nullptr)
    {

    }

    void activate(ECFExpressionVector &expression)
    {
        this->expression = &expression;
    }
};

struct ExternalEvaluator: SubKernel {
    const char* name() const { return evaluator->expression(); }

    ExternalEvaluator(Evaluator *evaluator)
    : evaluator(evaluator)
    {
        isconst = evaluator->isConst();
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
    }
    virtual ~ExternalEvaluator() {}

    Evaluator *evaluator;

    void setTime(double time, int t)
    {
        evaluator->getTime(t) = time;
    }

    void setFrequency(double frequency, int t)
    {
        evaluator->getFrequency(t) = frequency;
    }
};

template <size_t ndim> struct SetExpression;

template <> struct SetExpression<2> {
    template <typename Element>
    static void coordinates_node(Element &element, size_t n, size_t s, double &coordinateX, double &coordinateY, double &coordinateZ)
    {
        coordinateX = element.coords.node[n][0][s];
        coordinateY = element.coords.node[n][1][s];
    }

    template <typename Element>
    static void coordinates_gp(Element &element, size_t s, double &coordinateX, double &coordinateY, double &coordinateZ)
    {
        coordinateX = element.coords.gp[0][s];
        coordinateY = element.coords.gp[1][s];
    }
};

template <> struct SetExpression<3> {
    template <typename Element>
    static void coordinates_node(Element &element, size_t n, size_t s, double &coordinateX, double &coordinateY, double &coordinateZ)
    {
        coordinateX = element.coords.node[n][0][s];
        coordinateY = element.coords.node[n][1][s];
        coordinateZ = element.coords.node[n][2][s];
    }

    template <typename Element>
    static void coordinates_gp(Element &element, size_t s, double &coordinateX, double &coordinateY, double &coordinateZ)
    {
        coordinateX = element.coords.gp[0][s];
        coordinateY = element.coords.gp[1][s];
        coordinateZ = element.coords.gp[2][s];
    }
};

template <size_t ndim, class Element>
struct ExternalGPsExpression: ExternalEvaluator, Element {

    std::function<void(Element&, size_t&, size_t&, double)> setter;
    int t = omp_get_thread_num();
    double &temperature;
    double &coordinateX, &coordinateY, &coordinateZ;

    ExternalGPsExpression(Evaluator *evaluator, const std::function<void(Element&, size_t&, size_t&, double)> &setter)
    : ExternalEvaluator(evaluator), setter(setter),
      t(omp_get_thread_num()),
      temperature(evaluator->getTemperature(t)),
      coordinateX(evaluator->getCoordinateX(t)), coordinateY(evaluator->getCoordinateY(t)), coordinateZ(evaluator->getCoordinateZ(t))
    {

    }

    void simd(Element &element, size_t gp)
    {
        for (size_t s = 0; s < SIMD::size; ++s) {
            temperature = element.temperature.gp[s];
            SetExpression<ndim>::coordinates_gp(element, s, coordinateX, coordinateY, coordinateZ);
            setter(element, gp, s, this->evaluator->evaluate());
        }
    }
};

template <size_t ndim, class Element>
struct ExternalNodeExpression: ExternalEvaluator, Element {

    std::function<void(Element&, size_t&, size_t&, double)> setter;
    int t = omp_get_thread_num();
    double &temperature;
    double &coordinateX, &coordinateY, &coordinateZ;

    ExternalNodeExpression(Evaluator *evaluator, const std::function<void(Element&, size_t&, size_t&, double)> &setter)
    : ExternalEvaluator(evaluator), setter(setter),
      t(omp_get_thread_num()),
      temperature(evaluator->getTemperature(t)),
      coordinateX(evaluator->getCoordinateX(t)), coordinateY(evaluator->getCoordinateY(t)), coordinateZ(evaluator->getCoordinateZ(t))
    {

    }

    void simd(Element &element, size_t n)
    {
        for (size_t s = 0; s < SIMD::size; ++s) {
            temperature = element.temperature.node[n][s];
            SetExpression<ndim>::coordinates_node(element, n, s, coordinateX, coordinateY, coordinateZ);
            setter(element, n, s, this->evaluator->evaluate());
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_EXPRESSION_H_ */
