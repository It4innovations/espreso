
#ifndef SRC_CONFIG_ECF_PHYSICS_HEATTRANSFER_H_
#define SRC_CONFIG_ECF_PHYSICS_HEATTRANSFER_H_

#include "physics.h"
#include "physicssolver/loadstep.h"

namespace espreso {

struct ECF;

struct ConvectionConfiguration: public ECFDescription {

    enum class TYPE {
        USER,
        EXTERNAL_NATURAL,
        INTERNAL_NATURAL,
        EXTERNAL_FORCED,
        INTERNAL_FORCED
    };

    enum class VARIANT {
        VERTICAL_WALL,
        INCLINED_WALL,
        HORIZONTAL_CYLINDER,
        SPHERE,
        HORIZONTAL_PLATE_UP,
        HORIZONTAL_PLATE_DOWN,
        AVERAGE_PLATE,
        PARALLEL_PLATES,
        CIRCULAR_TUBE,
        TUBE,
        QUENCH,
        QUENCH_PARALLEL
    };

    enum class FLUID {
        AIR,
        WATER,
        ENGINE_OIL,
        TRANSFORMER_OIL,
        STEAM
    };

    TYPE type;
    VARIANT variant;
    FLUID fluid;

    ECFExpression heat_transfer_coefficient, external_temperature;
    ECFExpression wall_height, tilt_angle, diameter, plate_length, fluid_velocity, plate_distance, length, experimental_constant, volume_fraction, absolute_pressure;

    ConvectionConfiguration();
};

struct RadiationConfiguration: public ECFDescription {

    ECFExpression emissivity, external_temperature;

    RadiationConfiguration();
};

struct BioHeatSourceConfiguration: public ECFDescription {

    ECFExpression arteriar_blood_temperature;
    ECFExpression blood_specific_heat, blood_density;
    ECFExpression metabolic_heat_source;
    ECFExpression blood_perfusion;
    ECFExpression reference_temperature;
    ECFExpression physical_activity_scatter_factor;
    ECFExpression mu;

    BioHeatSourceConfiguration();
};

struct HeatTransferGlobalSettings {

    enum class STABILIZATION {
        SUPG = 0,
        CAU  = 1
    };

    enum class KERNEL {
        OLD,
        OPT,
        VEC
    };

    STABILIZATION stabilization;
    KERNEL kernel;
    double sigma;
    bool init_temp_respect_bc, diffusion_split;

    HeatTransferGlobalSettings(ECFObject *ecfdescription);
};

struct HumanThermoregulationSystem: public ECFDescription {

    enum class ACTIVITY_LEVEL_UNIT {
        STANDING        = 0,
        EASY_WORK       = 1,
        HARD_WORK       = 2,
        WALKING_0_PRCT  = 3,
        WALKING_5_PRCT  = 4,
        WALKING_15_PRCT = 5,
        TEACHER         = 6
    };

    ACTIVITY_LEVEL_UNIT activity_level_unit;
    HumanThermoregulationSystem();
};

struct HeatTransferLoadStepConfiguration: public HeatTransferLoadStepSolverConfiguration {

    bool update_initial_temperature;

    std::map<std::string, ECFExpression> temperature, heat_source, heat_flux, heat_flow;
    std::map<std::string, ECFExpressionVector> translation_motions;
    std::map<std::string, BioHeatSourceConfiguration> bio_heat;
    std::map<std::string, ConvectionConfiguration> convection;
    std::map<std::string, RadiationConfiguration> diffuse_radiation;
    HumanThermoregulationSystem human_thermoregulation_system;

    HeatTransferLoadStepConfiguration();
};

struct HeatTransferOutputSettings: public virtual ECFDescription {

    static void addMonitorableProperties(ECFMetaData &metadata, const ECF *root);

    bool temperature, translation_motions, gradient, flux, phase, latent_heat, htc;

    void basic() {
        temperature = translation_motions = true;
        gradient = flux = phase = latent_heat = htc = false;
    }
    void all() {
        temperature = translation_motions = gradient = flux = phase = latent_heat = htc = true;
    }

    static void activate();

    HeatTransferOutputSettings();

protected:
    static bool _activated;
};

struct HeatTransferConfiguration: public PhysicsConfiguration, public HeatTransferGlobalSettings {

    std::map<size_t, HeatTransferLoadStepConfiguration> load_steps_settings;

    HeatTransferConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_HEATTRANSFER_H_ */
