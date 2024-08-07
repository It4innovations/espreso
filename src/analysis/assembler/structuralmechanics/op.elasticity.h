
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELASTICITY_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/op.coordinatesystem.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct Elasticity: SubKernel {
    const char* name() const { return "Elasticity"; }

    Elasticity()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN),
      configuration(nullptr),
      coordinateSystem(nullptr),
      type(CoordinateSystemConfiguration::TYPE::CARTESIAN),
      rotated(false),
      constRotation(true)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
    }

    void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour, const LinearElasticPropertiesConfiguration *elasticity, const CoordinateSystemConfiguration *coordinateSystem)
    {
        this->behaviour = behaviour;
        this->configuration = elasticity;
        this->coordinateSystem = coordinateSystem;
        this->isconst = !(elasticity->needCoordinates() || elasticity->needTemperature()) && coordinateSystem->isConst();
        this->type = coordinateSystem->type;
        this->rotated = coordinateSystem->isRotated() && elasticity->model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC;
        this->constRotation = coordinateSystem->type == CoordinateSystemConfiguration::TYPE::CARTESIAN && coordinateSystem->isConst();
        this->isactive = !this->isconst || !this->constRotation;
    }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
    const LinearElasticPropertiesConfiguration *configuration;
    const CoordinateSystemConfiguration *coordinateSystem;
    CoordinateSystemConfiguration::TYPE type;
    bool rotated, constRotation;
};

template <size_t ndim> struct ElasticityKernel;

template <> struct ElasticityKernel<2>: Elasticity {
    ElasticityKernel(const Elasticity &base): Elasticity(base) {}

    template <typename Element>
    void simd(Element &element)
    {
//        CoordinateSystem<2, 2>::simd(element, type, false);
        simd(element, 0);
    }


    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (configuration->model) {
        case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
            switch (behaviour) {
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN: {
                SIMD ex = element.ecf.youngModulus[0];
                SIMD mi = element.ecf.poissonRatio[0];
                SIMD C1 = load1(1.);
                SIMD C2 = load1(2.);
                SIMD k = ex * (C1 - mi) / ((C1 + mi) * (C1 - C2 * mi));
                if (rotated) {
                    element.ecf.elasticity[0] = element.ecf.elasticity[4] = k;
                    element.ecf.elasticity[1] = k * (mi / (C1 - mi));
                    element.ecf.elasticity[8] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
                } else {
                    element.elasticity[0] = element.elasticity[4] = k;
                    element.elasticity[1] = k * (mi / (C1 - mi));
                    element.elasticity[8] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
                }
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS: {
                SIMD ex = element.ecf.youngModulus[0];
                SIMD mi = element.ecf.poissonRatio[0];
                SIMD C1 = load1(1.);
                SIMD C2 = load1(.5);
                SIMD k = ex / (C1 - mi * mi);
                if (rotated) {
                    element.ecf.elasticity[0] = element.ecf.elasticity[4] = k;
                    element.ecf.elasticity[1] = k * mi;
                    element.ecf.elasticity[8] = k * ((C1 -  mi) * C2);
                } else {
                    element.elasticity[0] = element.elasticity[4] = k;
                    element.elasticity[1] = k * mi;
                    element.elasticity[8] = k * ((C1 -  mi) * C2);
                }
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC: {
                SIMD ex = element.ecf.youngModulus[0];
                SIMD mi = element.ecf.poissonRatio[0];
                SIMD C05 = load1(.5);
                SIMD C1 = load1(1.);
                SIMD C2 = load1(2.);
                SIMD k = ex / ((C1 + mi) * (C1 - C2 * mi));
                if (rotated) {
                    element.ecf.elasticity[0] = element.ecf.elasticity[5] = element.ecf.elasticity[10] = k * (C1 - mi);
                    element.ecf.elasticity[1] = element.ecf.elasticity[2] = element.ecf.elasticity[6] = k * mi;
                    element.ecf.elasticity[15] = k * (C1 - C2 * mi) * C05;
                } else {
                    element.elasticity[0] = element.elasticity[5] = element.elasticity[10] = k * (C1 - mi);
                    element.elasticity[1] = element.elasticity[2] = element.elasticity[6] = k * mi;
                    element.elasticity[15] = k * (C1 - C2 * mi) * C05;
                }
            } break;
            } break;
        case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
            switch (behaviour) {
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN: {
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS: {
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC: {
            } break;
            } break;
        case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
            switch (behaviour) {
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN: {
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS: {
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC: {
            } break;
            } break;
        }

        if (rotated) {
            // make the elasticity symmetric
            element.ecf.elasticity[3] = element.ecf.elasticity[1];
            element.ecf.elasticity[6] = element.ecf.elasticity[2]; element.ecf.elasticity[7] = element.ecf.elasticity[5];
        } else {
            // make the elasticity symmetric
            element.elasticity[3] = element.elasticity[1];
            element.elasticity[6] = element.elasticity[2]; element.elasticity[7] = element.elasticity[5];
        }
    }
};

template <> struct ElasticityKernel<3>: Elasticity {
    ElasticityKernel(const Elasticity &base): Elasticity(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        CoordinateSystem<3, 2>::simd(element, type, false);
        simd(element, 0);
    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (configuration->model) {
        case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
        {
            // C
            // 0 1 1 _ _ _
            //   0 1 _ _ _
            //     0 _ _ _
            //       2 _ _
            //         2 _
            //           2
            SIMD ex = element.ecf.youngModulus[0];
            SIMD mi = element.ecf.poissonRatio[0];
            SIMD C05 = load1(.5);
            SIMD C1 = load1(1.);
            SIMD C2 = load1(2.);
            SIMD ee = ex / ((C1 + mi) * (C1 - C2 * mi));
//            if (rotated) {
//                element.ecf.elasticity[0]  = element.ecf.elasticity[7]  = element.ecf.elasticity[14] = ee * (C1 - mi);
//                element.ecf.elasticity[1]  = element.ecf.elasticity[2]  = element.ecf.elasticity[8]  = ee * mi;
//                element.ecf.elasticity[21] = element.ecf.elasticity[28] = element.ecf.elasticity[35] = ee * (C05 - mi);
//            } else {
                element.elasticity[ 0] = element.elasticity[ 7] = element.elasticity[14] = ee * (C1 - mi);
                element.elasticity[ 1] = element.elasticity[ 2] = element.elasticity[ 8]  = ee * mi;
                element.elasticity[21] = element.elasticity[28] = element.elasticity[35] = ee * (C05 - mi);
//            }
        } break;
        case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
        {
            // C
            // 0 1 2 _ _ _
            //   3 4 _ _ _
            //     5 _ _ _
            //       6 _ _
            //         7 _
            //           8
            SIMD ex = element.ecf.youngModulus[0];
            SIMD ey = element.ecf.youngModulus[1];
            SIMD ez = element.ecf.youngModulus[2];
            SIMD miXY = element.ecf.poissonRatio[0];
            SIMD miXZ = element.ecf.poissonRatio[1];
            SIMD miYZ = element.ecf.poissonRatio[2];

            SIMD miYX = miXY * ey / ex;
            SIMD miZY = miYZ * ez / ey;
            SIMD miZX = miXZ * ex / ez;

            SIMD C1 = load1(1);
            SIMD ksi = C1 / (C1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ));

            if (rotated) {
                element.ecf.elasticity[0 ] = ksi * ex * (C1   - miYZ * miZY);
                element.ecf.elasticity[1 ] = ksi * ey * (miXY + miXZ * miZY);
                element.ecf.elasticity[2 ] = ksi * ez * (miXZ + miYZ * miXY);
                element.ecf.elasticity[7 ] = ksi * ey * (C1   - miXZ * miZX);
                element.ecf.elasticity[8 ] = ksi * ez * (miYZ + miYX * miXZ);
                element.ecf.elasticity[14] = ksi * ez * (C1   - miYX * miXY);
                element.ecf.elasticity[21] = element.ecf.shearModulus[0];
                element.ecf.elasticity[28] = element.ecf.shearModulus[1];
                element.ecf.elasticity[35] = element.ecf.shearModulus[2];

                SIMD C05 = load1(0.5);
                CoordinateSystem<3, 2>::simd(element, type, constRotation);
                switch (type) { // always rotated coordinate system
                case CoordinateSystemConfiguration::TYPE::CARTESIAN: {
                    SIMD cosx  = element.rotation.cossin[ 0];
                    SIMD cosy  = element.rotation.cossin[ 1];
                    SIMD cosz  = element.rotation.cossin[ 2];
                    SIMD sinx  = element.rotation.cossin[ 3];
                    SIMD siny  = element.rotation.cossin[ 4];
                    SIMD sinz  = element.rotation.cossin[ 5];
                    SIMD cos2x = element.rotation.cossin[ 6];
                    SIMD cos2y = element.rotation.cossin[ 7];
                    SIMD cos2z = element.rotation.cossin[ 8];
                    SIMD sin2x = element.rotation.cossin[ 9];
                    SIMD sin2y = element.rotation.cossin[10];
                    SIMD sin2z = element.rotation.cossin[11];

                    SIMD T00 = cosx * cosx * cosz * cosz - sin2x * sin2z * cosy * C05 + cosy * cosy * sinx * sinx * sinz * sinz;
                    SIMD T01 = cosz * cosz * sinx * sinx + sin2x * sin2z * cosy * C05 + cosx * cosx * cosy * cosy * sinz * sinz;
                    SIMD T02 = siny * siny * sinz * sinz;
                    SIMD T03 = sin2z * sinx * siny + sin2y * cosx * sinz * sinz;
                    SIMD T04 = sin2z * cosx * siny - sin2y * sinx * sinz * sinz;
                    SIMD T05 = sin2x * cosz * cosz + cos2x * sin2z * cosy - sin2x * cosy * cosy * sinz * sinz;
                    SIMD T10 = cosx * cosx * sinz * sinz + sin2x * sin2z * cosy * C05 + cosy * cosy * cosz * cosz * sinx * sinx;
                    SIMD T11 = sinx * sinx * sinz * sinz - sin2x * sin2z * cosy * C05 + cosx * cosx * cosy * cosy * cosz * cosz;
                    SIMD T12 = cosz * cosz * siny * siny;
                    SIMD T13 = -sin2z * sinx * siny + sin2y * cosx * cosz * cosz;
                    SIMD T14 = -sin2z * cosx * siny - sin2y * cosz * cosz * sinx;
                    SIMD T15 = sin2x * sinz * sinz - cos2x * sin2z * cosy - sin2x * cosy * cosy * cosz * cosz;
                    SIMD T20 = sinx * sinx * siny * siny;
                    SIMD T21 = cosx * cosx * siny * siny;
                    SIMD T22 = cosy * cosy;
                    SIMD T23 = -sin2y * cosx;
                    SIMD T24 = sin2y * sinx;
                    SIMD T25 = -sin2x * siny * siny;
                    SIMD T30 = -C05 * sin2x * siny * sinz - sin2y * cosz * sinx * sinx * C05;
                    SIMD T31 = sin2x * siny * sinz * C05 - sin2y * cosx * cosx * cosz * C05;
                    SIMD T32 = sin2y * cosz * C05;
                    SIMD T33 = cos2y * cosx * cosz - cosy * sinx * sinz;
                    SIMD T34 = -cos2y * cosz * sinx - cosx * cosy * sinz;
                    SIMD T35 = cos2x * siny * sinz + sin2x * sin2y * cosz * C05;
                    SIMD T40 = sin2x * cosz * siny * C05 - sin2y * sinx * sinx * sinz * C05;
                    SIMD T41 = -C05 * sin2x * cosz * siny - sin2y * cosx * cosx * sinz * C05;
                    SIMD T42 = sin2y * sinz * C05;
                    SIMD T43 = cos2y * cosx * sinz + cosy * cosz * sinx;
                    SIMD T44 = -cos2y * sinx * sinz + cosx * cosy * cosz;
                    SIMD T45 = -cos2x * cosz * siny + sin2x * sin2y * sinz * C05;
                    SIMD T50 = -C05 * sin2z * cosx * cosx - cos2z * sin2x * cosy * C05 + sin2z * cosy * cosy * sinx * sinx * C05;
                    SIMD T51 = -C05 * sin2z * sinx * sinx + cos2z * sin2x * cosy * C05 + sin2z * cosx * cosx * cosy * cosy * C05;
                    SIMD T52 = sin2z * siny * siny * C05;
                    SIMD T53 = cos2z * sinx * siny + sin2y * sin2z * cosx * C05;
                    SIMD T54 = cos2z * cosx * siny - sin2y * sin2z * sinx * C05;
                    SIMD T55 = -C05 * sin2x * sin2z + cos2x * cos2z * cosy - sin2x * sin2z * cosy * cosy * C05;

                    SIMD e00 = element.ecf.elasticity[0];
                    SIMD e01 = element.ecf.elasticity[1], e11 = element.ecf.elasticity[7];
                    SIMD e02 = element.ecf.elasticity[2], e12 = element.ecf.elasticity[8], e22 = element.ecf.elasticity[14];
                    SIMD e33 = element.ecf.elasticity[21];
                    SIMD e44 = element.ecf.elasticity[28];
                    SIMD e55 = element.ecf.elasticity[35];

                    SIMD a = T00 * e00 + T01 * e01 + T02 * e02;
                    SIMD b = T00 * e01 + T01 * e11 + T02 * e12;
                    SIMD c = T00 * e02 + T01 * e12 + T02 * e22;
                    SIMD d = T03 * e33;
                    SIMD e = T04 * e44;
                    SIMD f = T05 * e55;
                    element.elasticity[ 0] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
                    element.elasticity[ 1] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
                    element.elasticity[ 2] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
                    element.elasticity[ 3] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
                    element.elasticity[ 4] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[ 5] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T10 * e00 + T11 * e01 + T12 * e02;
                    b = T10 * e01 + T11 * e11 + T12 * e12;
                    c = T10 * e02 + T11 * e12 + T12 * e22;
                    d = T13 * e33;
                    e = T14 * e44;
                    f = T15 * e55;
                    element.elasticity[ 7] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
                    element.elasticity[ 8] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
                    element.elasticity[ 9] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
                    element.elasticity[10] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[11] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T20 * e00 + T21 * e01 + T22 * e02;
                    b = T20 * e01 + T21 * e11 + T22 * e12;
                    c = T20 * e02 + T21 * e12 + T22 * e22;
                    d = T23 * e33;
                    e = T24 * e44;
                    f = T25 * e55;
                    element.elasticity[14] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
                    element.elasticity[15] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
                    element.elasticity[16] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[17] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T30 * e00 + T31 * e01 + T32 * e02;
                    b = T30 * e01 + T31 * e11 + T32 * e12;
                    c = T30 * e02 + T31 * e12 + T32 * e22;
                    d = T33 * e33;
                    e = T34 * e44;
                    f = T35 * e55;
                    element.elasticity[21] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
                    element.elasticity[22] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[23] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T40 * e00 + T41 * e01 + T42 * e02;
                    b = T40 * e01 + T41 * e11 + T42 * e12;
                    c = T40 * e02 + T41 * e12 + T42 * e22;
                    d = T43 * e33;
                    e = T44 * e44;
                    f = T45 * e55;
                    element.elasticity[28] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[29] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T50 * e00 + T51 * e01 + T52 * e02;
                    b = T50 * e01 + T51 * e11 + T52 * e12;
                    c = T50 * e02 + T51 * e12 + T52 * e22;
                    d = T53 * e33;
                    e = T54 * e44;
                    f = T55 * e55;
                    element.elasticity[35] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;
                } break;
                case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {
                    SIMD cosz  = element.rotation.cossin[ 2];
                    SIMD sinz  = element.rotation.cossin[ 5];
                    SIMD cos2z = element.rotation.cossin[ 8];
                    SIMD sin2z = element.rotation.cossin[11];

                    SIMD T00 =  cosz * cosz;
                    SIMD T01 =  sinz * sinz;
                    SIMD T05 =  sin2z;
                    SIMD T10 =  sinz * sinz;
                    SIMD T11 =  cosz * cosz;
                    SIMD T15 = -sin2z;
                    SIMD T22 =  load1(1);
                    SIMD T33 =  cosz;
                    SIMD T34 = -sinz;
                    SIMD T43 =  sinz;
                    SIMD T44 =  cosz;
                    SIMD T50 = -C05 * sin2z;
                    SIMD T51 =  sin2z * C05;
                    SIMD T55 =  cos2z;

                    SIMD e00 = element.ecf.elasticity[0];
                    SIMD e01 = element.ecf.elasticity[1], e11 = element.ecf.elasticity[7];
                    SIMD e02 = element.ecf.elasticity[2], e12 = element.ecf.elasticity[8], e22 = element.ecf.elasticity[14];
                    SIMD e33 = element.ecf.elasticity[21];
                    SIMD e44 = element.ecf.elasticity[28];
                    SIMD e55 = element.ecf.elasticity[35];

                    SIMD a = T00 * e00 + T01 * e01;
                    SIMD b = T00 * e01 + T01 * e11;
                    SIMD c = T00 * e02 + T01 * e12;
                    SIMD d;
                    SIMD e;
                    SIMD f = T05 * e55;
                    element.elasticity[ 0] = a * T00 + b * T01                               + f * T05;
                    element.elasticity[ 1] = a * T10 + b * T11                               + f * T15;
                    element.elasticity[ 2] =                     c * T22                              ;
                    element.elasticity[ 3] = load1(0)                                                 ;
                    element.elasticity[ 4] = load1(0)                                                 ;
                    element.elasticity[ 5] = a * T50 + b * T51                               + f * T55;

                    a = T10 * e00 + T11 * e01;
                    b = T10 * e01 + T11 * e11;
                    c = T10 * e02 + T11 * e12;
                    f = T15 * e55;
                    element.elasticity[ 7] = a * T10 + b * T11 +                             + f * T15;
                    element.elasticity[ 8] =                     c * T22                              ;
                    element.elasticity[ 9] = load1(0)                                                 ;
                    element.elasticity[10] = load1(0)                                                 ;
                    element.elasticity[11] = a * T50 + b * T51                               + f * T55;

                    a = T22 * e02;
                    b = T22 * e12;
                    c = T22 * e22;
                    element.elasticity[14] =                     c * T22                              ;
                    element.elasticity[15] = load1(0)                                                 ;
                    element.elasticity[16] = load1(0)                                                 ;
                    element.elasticity[17] = a * T50 + b * T51                                        ;

                    d = T33 * e33;
                    e = T34 * e44;
                    element.elasticity[21] =                               d * T33 + e * T34          ;
                    element.elasticity[22] =                               d * T43 + e * T44          ;
                    element.elasticity[23] = load1(0)                                                 ;

                    d = T43 * e33;
                    e = T44 * e44;
                    element.elasticity[28] =                               d * T43 + e * T44          ;
                    element.elasticity[29] = load1(0)                                                 ;

                    a = T50 * e00 + T51 * e01;
                    b = T50 * e01 + T51 * e11;
                    f = T55 * e55;
                    element.elasticity[35] = a * T50 + b * T51                               + f * T55;
                } break;
                case CoordinateSystemConfiguration::TYPE::SPHERICAL: {
                    SIMD cosx  = element.rotation.cossin[ 0];
                    SIMD cosy  = element.rotation.cossin[ 1];
                    SIMD cosz  = element.rotation.cossin[ 2];
                    SIMD siny  = element.rotation.cossin[ 4];
                    SIMD sinz  = element.rotation.cossin[ 5];
                    SIMD cos2x = element.rotation.cossin[ 6];
                    SIMD cos2y = element.rotation.cossin[ 7];
                    SIMD cos2z = element.rotation.cossin[ 8];
                    SIMD sin2y = element.rotation.cossin[10];
                    SIMD sin2z = element.rotation.cossin[11];

                    SIMD T00 =  cosx * cosx * cosz * cosz;
                    SIMD T01 =  cosx * cosx * cosy * cosy * sinz * sinz;
                    SIMD T02 =  siny * siny * sinz * sinz;
                    SIMD T03 =  sin2y * cosx * sinz * sinz;
                    SIMD T04 =  sin2z * cosx * siny;
                    SIMD T05 =  cos2x * sin2z * cosy;
                    SIMD T10 =  cosx * cosx * sinz * sinz;
                    SIMD T11 =  cosx * cosx * cosy * cosy * cosz * cosz;
                    SIMD T12 =  cosz * cosz * siny * siny;
                    SIMD T13 =  sin2y * cosx * cosz * cosz;
                    SIMD T14 = -sin2z * cosx * siny;
                    SIMD T15 = -cos2x * sin2z * cosy;
                    SIMD T21 =  cosx * cosx * siny * siny;
                    SIMD T22 =  cosy * cosy;
                    SIMD T23 = -sin2y * cosx;
                    SIMD T31 = -sin2y * cosx * cosx * cosz * C05;
                    SIMD T32 =  sin2y * cosz * C05;
                    SIMD T33 =  cos2y * cosx * cosz;
                    SIMD T34 = -cosx * cosy * sinz;
                    SIMD T35 =  cos2x * siny * sinz;
                    SIMD T41 = -sin2y * cosx * cosx * sinz * C05;
                    SIMD T42 =  sin2y * sinz * C05;
                    SIMD T43 =  cos2y * cosx * sinz;
                    SIMD T44 =  cosx * cosy * cosz;
                    SIMD T45 = -cos2x * cosz * siny;
                    SIMD T50 = -C05 * sin2z * cosx * cosx;
                    SIMD T51 =  sin2z * cosx * cosx * cosy * cosy * C05;
                    SIMD T52 =  sin2z * siny * siny * C05;
                    SIMD T53 =  sin2y * sin2z * cosx * C05;
                    SIMD T54 =  cos2z * cosx * siny;
                    SIMD T55 =  cos2x * cos2z * cosy;

                    SIMD e00 = element.ecf.elasticity[0];
                    SIMD e01 = element.ecf.elasticity[1], e11 = element.ecf.elasticity[7];
                    SIMD e02 = element.ecf.elasticity[2], e12 = element.ecf.elasticity[8], e22 = element.ecf.elasticity[14];
                    SIMD e33 = element.ecf.elasticity[21];
                    SIMD e44 = element.ecf.elasticity[28];
                    SIMD e55 = element.ecf.elasticity[35];

                    SIMD a = T00 * e00 + T01 * e01 + T02 * e02;
                    SIMD b = T00 * e01 + T01 * e11 + T02 * e12;
                    SIMD c = T00 * e02 + T01 * e12 + T02 * e22;
                    SIMD d = T03 * e33;
                    SIMD e = T04 * e44;
                    SIMD f = T05 * e55;
                    element.elasticity[ 0] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
                    element.elasticity[ 1] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
                    element.elasticity[ 2] =           b * T21 + c * T22 + d * T23                    ;
                    element.elasticity[ 3] =           b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
                    element.elasticity[ 4] =           b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[ 5] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T10 * e00 + T11 * e01 + T12 * e02;
                    b = T10 * e01 + T11 * e11 + T12 * e12;
                    c = T10 * e02 + T11 * e12 + T12 * e22;
                    d = T13 * e33;
                    e = T14 * e44;
                    f = T15 * e55;
                    element.elasticity[ 7] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
                    element.elasticity[ 8] =           b * T21 + c * T22 + d * T23                    ;
                    element.elasticity[ 9] =           b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
                    element.elasticity[10] =           b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[11] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T21 * e01 + T22 * e02;
                    b = T21 * e11 + T22 * e12;
                    c = T21 * e12 + T22 * e22;
                    d = T23 * e33;
                    element.elasticity[14] =           b * T21 + c * T22 + d * T23                    ;
                    element.elasticity[15] =           b * T31 + c * T32 + d * T33                    ;
                    element.elasticity[16] =           b * T41 + c * T42 + d * T43                    ;
                    element.elasticity[17] = a * T50 + b * T51 + c * T52 + d * T53                    ;

                    a = T31 * e01 + T32 * e02;
                    b = T31 * e11 + T32 * e12;
                    c = T31 * e12 + T32 * e22;
                    d = T33 * e33;
                    e = T34 * e44;
                    f = T35 * e55;
                    element.elasticity[21] =           b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
                    element.elasticity[22] =           b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[23] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T41 * e01 + T42 * e02;
                    b = T41 * e11 + T42 * e12;
                    c = T41 * e12 + T42 * e22;
                    d = T43 * e33;
                    e = T44 * e44;
                    f = T45 * e55;
                    element.elasticity[28] =           b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
                    element.elasticity[29] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

                    a = T50 * e00 + T51 * e01 + T52 * e02;
                    b = T50 * e01 + T51 * e11 + T52 * e12;
                    c = T50 * e02 + T51 * e12 + T52 * e22;
                    d = T53 * e33;
                    e = T54 * e44;
                    f = T55 * e55;
                    element.elasticity[35] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;
                } break;
                }
            } else {
                element.elasticity[0 ] = ksi * ex * (C1   - miYZ * miZY);
                element.elasticity[1 ] = ksi * ey * (miXY + miXZ * miZY);
                element.elasticity[2 ] = ksi * ez * (miXZ + miYZ * miXY);
                element.elasticity[7 ] = ksi * ey * (C1   - miXZ * miZX);
                element.elasticity[8 ] = ksi * ez * (miYZ + miYX * miXZ);
                element.elasticity[14] = ksi * ez * (C1   - miYX * miXY);
                element.elasticity[21] = element.ecf.shearModulus[0];
                element.elasticity[28] = element.ecf.shearModulus[1];
                element.elasticity[35] = element.ecf.shearModulus[2];
            }
        } break;
        case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
        {

        } break;
        }

        // make the elasticity symmetric
        element.elasticity[ 6] = element.elasticity[ 1];
        element.elasticity[12] = element.elasticity[ 2]; element.elasticity[13] = element.elasticity[ 8];
        element.elasticity[18] = element.elasticity[ 3]; element.elasticity[19] = element.elasticity[ 9]; element.elasticity[20] = element.elasticity[15];
        element.elasticity[24] = element.elasticity[ 4]; element.elasticity[25] = element.elasticity[10]; element.elasticity[26] = element.elasticity[16]; element.elasticity[27] = element.elasticity[22];
        element.elasticity[30] = element.elasticity[ 5]; element.elasticity[31] = element.elasticity[11]; element.elasticity[32] = element.elasticity[17]; element.elasticity[33] = element.elasticity[23]; element.elasticity[34] = element.elasticity[29];
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_ELASTICITY_H_ */
