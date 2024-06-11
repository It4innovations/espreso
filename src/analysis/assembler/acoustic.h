
#ifndef SRC_ANALYSIS_ASSEMBLER_ACOUSTIC_H_
#define SRC_ANALYSIS_ASSEMBLER_ACOUSTIC_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"
#include "assembler.h"
#include "config/ecf/physics/acoustic.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"

#include <vector>

namespace espreso {

struct AcousticConfiguration;
struct AcousticLoadStepConfiguration;
struct Harmonic;

class Acoustic: public Assembler
{
public:
    struct NGP {
        static const size_t POINT1 = 1;

        static const size_t LINE2 = 2;
        static const size_t LINE3 = 3;

        static const size_t TRIANGLE3 = 6;
        static const size_t TRIANGLE6 = 6;
        static const size_t SQUARE4   = 4;
        static const size_t SQUARE8   = 9;

        static const size_t TETRA4    = 4;
        static const size_t TETRA10   = 15;
        static const size_t PYRAMID5  = 8;
        static const size_t PYRAMID13 = 14;
        static const size_t PRISMA6   = 9;
        static const size_t PRISMA15  = 9;
        static const size_t HEXA8     = 8;
        static const size_t HEXA20    = 8;
    };

    Acoustic(Acoustic *previous, AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration);

    void analyze(const step::Step &step);

    void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *dirichlet);
    void evaluate(const step::Step &step, step::Frequency &frequency, Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *dirichlet);
    void updateSolution(const step::Step &step, Vector_Base<double> *rex, Vector_Base<double> *imx);

    AcousticConfiguration &settings;
    AcousticLoadStepConfiguration &configuration;

//    struct ParametersAcousticPressure {
//        ElementParameter<enodes> node;
//        ElementGPsExternalParameter<egps> gp;
//
//        struct {
//            ElementGPsExternalParameter<enodes> node;
//            ElementGPsExternalParameter<egps> gp;
//        } initial;
//    } acoustic_pressure;
//
//    struct {
//        ElementParameter<egps> weight;
//        ElementParameter<enodes * egps> N;
//        ElementParameter<edim * enodes * egps> dN;
//
//        ElementParameter<egps> jacobiDeterminant;
//        ElementParameter<ndim * ndim * egps> jacobiInversion;
//        ElementParameter<edim * enodes * egps> dND;
//
//        struct {
//            BoundaryParameter<egps> weight;
//            BoundaryParameter<enodes * egps> N;
//            BoundaryParameter<edim * enodes * egps> dN;
//
//            BoundaryParameter<egps> jacobian;
//        } boundary;
//    } integration;
//
//    struct {
//        ElementParameter<ndim * enodes> node;
//        ElementParameter<ndim * egps> gp;
//        struct {
//            BoundaryParameter<ndim * enodes> node;
//            BoundaryParameter<ndim * egps> gp;
//        } boundary;
//    } coords;

//    struct {
//        BoundaryExternalParameter<enodes> node;
//    } pressure;
//
//    struct {
//        BoundaryExternalParameter<egps> gp;
//    } normalAcceleration, impedance, q, proj_acceleration;
//
//    struct {
//        BoundaryExternalParameter<ndim * egps> gp;
//    } acceleration, normals;
//
//    struct {
//        ElementGPsExternalParameter<egps> density, speed_of_sound;
//    } material;
//
//    struct {
//        BoundaryExternalParameter<enodes> node;
//    } pointSource;

//    struct {
//        ElementGPsExternalParameter<egps> gp;
//    } monopoleSource;
//
//    struct {
//        ElementGPsExternalParameter<ndim * egps> gp;
//    } dipoleSource;

    struct Results {
        static NodeData *pressure, *initialPressure;
    };
protected:
    void elements(const step::Step &step, SubKernel::Action action, size_t interval) { }
    void boundary(const step::Step &step, SubKernel::Action action, size_t region, size_t interval) { }
    void nodes(const step::Step &step, SubKernel::Action action, size_t region, size_t interval) { }

    void initParameters();
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_ACOUSTIC_H_ */
