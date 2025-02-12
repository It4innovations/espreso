
#ifndef SRC_PHYSICS_KERNELS_MOVER_HEATTRANSFER_MOVER_H_
#define SRC_PHYSICS_KERNELS_MOVER_HEATTRANSFER_MOVER_H_

#include "mover.h"
#include "moverparameter.h"

namespace espreso {

class ConvectionConfiguration;
class RadiationConfiguration;
class PhysicsConfiguration;
class HeatTransferGlobalSettings;
class HeatTransferLoadStepConfiguration;
class HeatTransferConfiguration;
class BioHeatSourceConfiguration;

struct InputBioHeat {
    std::map<std::string, BioHeatSourceConfiguration> *ecf;
    Evaluator::Params *params;

    InputBioHeat(): ecf(NULL), params(NULL) {}

    bool isConst() const { return false; }
};

template <>
struct Move<InputBioHeat, ElementNodeValues>: public Mover {
    InputBioHeat from; ElementNodeValues bioHeat; bool moved;

    Move(const InputBioHeat &from, const ElementNodeValues &htc);

    void operator()();
    void now() { operator()(); }
};

struct InputBioHeatDerivation {
    std::map<std::string, BioHeatSourceConfiguration> *ecf;
    Evaluator::Params *params;

    InputBioHeatDerivation(): ecf(NULL), params(NULL) {}

    bool isConst() const { return false; }
};

template <>
struct Move<InputBioHeatDerivation, ElementNodeValues>: public Mover {
    InputBioHeatDerivation from; ElementNodeValues bioHeatDerivation; bool moved;

    Move(const InputBioHeatDerivation &from, const ElementNodeValues &htc);

    void operator()();
    void now() { operator()(); }
};

struct InputConvection {
    ConvectionConfiguration *ecf;
    Evaluator::Params *params;

    InputConvection(): ecf(NULL), params(NULL) {}

    bool isConst() const { return false; }
};

template <>
struct Move<InputConvection, ElementNodeValues>: public Mover {
    InputConvection from; ElementNodeValues htc; bool moved;

    Move(const InputConvection &from, const ElementNodeValues &htc);

    void operator()();
    void now() { operator()(); }
};

struct HeatTransferElementIterator: public ElementIterator {

    double activity_level_unit;
    MoverCoordinates coordinates;

    MoverFullParameter<InputExpressionMap      , OutputNodes>    temperature;
    MoverFullParameter<InputExpressionMap      , OutputElements> thickness;
    MoverFullParameter<InputExpressionVectorMap, OutputElements> motion;

    MoverInputParameter<InputExpressionMap>     initialTemperature;
    MoverInputParameter<InputExpressionMap>     heatSource;
    MoverInputParameter<InputBioHeat>           bioHeat;
    MoverInputParameter<InputBioHeatDerivation> bioHeatDerivation;

    MoverOutputParameter<OutputElements> phase;
    MoverOutputParameter<OutputElements> latentHeat;
    MoverOutputParameter<OutputElements> gradient;
    MoverOutputParameter<OutputElements> flux;
    MoverOutputParameter<OutputNodes> htc;

    HeatTransferElementIterator(HeatTransferElementIterator *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration, int dimension);
};

struct HeatTransferBoundaryIterator: public BoundaryIterator {
    double regionArea;
    bool radiation, convection;

    MoverCoordinates coordinates;
    MoverReducedInputParameter<OutputNodes> temperature;
    MoverReducedInputParameter<OutputElements> thickness;

    MoverInputParameter<InputExpression> emissivity;
    MoverInputParameter<InputExpression> extemperature;
    MoverInputParameter<InputExpression> heatflow;
    MoverInputParameter<InputExpression> heatflux;
    MoverBoundaryParameter<InputConvection, OutputNodes> htc;

    HeatTransferBoundaryIterator(BoundaryRegionStore *region, HeatTransferElementIterator &iterator, HeatTransferLoadStepConfiguration &configuration, int dimension);

    bool hasSettings();
};

}


#endif /* SRC_PHYSICS_KERNELS_MOVER_HEATTRANSFER_MOVER_H_ */
