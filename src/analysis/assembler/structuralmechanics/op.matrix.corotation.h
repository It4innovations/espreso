
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_COROTATION_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_COROTATION_H_

#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct MatrixCorotation: SubKernel {
    const char* name() const { return "MatrixCorotation"; }

    MatrixCorotation()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate()
    {
        this->isactive = 1;
    }
};

template <size_t nodes, size_t ndim> struct MatrixCorotationKernel;

template <size_t nodes>
struct MatrixCorotationKernel<nodes, 2>: MatrixCorotation {
    MatrixCorotationKernel(const MatrixCorotation &base): MatrixCorotation(base) {}

    template <typename Element>
    void simd(Element &element)
    {

    }
};

template <size_t nodes>
struct MatrixCorotationKernel<nodes, 3>: MatrixCorotation {
    MatrixCorotationKernel(const MatrixCorotation &base): MatrixCorotation(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        // x = coo + disp

        // teziste v referencnim
        // dXi = 3 x nodes : derivace bazovych funkci podle Xi

        // J = x0 * dXi  // kontrola na objem
        // dND = invJ * dXi

        // S = dND[i]
        // [  0 -3  2 ]
        // [  3  0 -1 ]
        // [ -2  1  0 ]

        // A_transp

        // D = disp * dND'
        // F = eye(3) + D

        // R = F * (F' * F)^(-1/2)
        // Re = blkdial(R,...)

        // teziste v coords + disp
        // teziste v pocatecnich

        // pro kazdy uzel vzdalenost od teziste v coords + disp
        // pro kazdy uzel vzdalenost od teziste v pocatecnich

        // S_bar = -R' * vektor od teziste (deformovane) * R

        // u_lie = Re' * vektor od teziste (deformovane) - vektor od teziste (nedef.)
        // f_ile = Ke * u_lie

        // internal force computation
        // G_bar = inv(A_trans * S_bar) * A_transp
        // P_bar = eye(24) - S_bar * G_bar

        // f -= fie
        // fie = Re * P_bar' * f_ile

        // tangent stiffness matrix
        // fiep = P_bar' * f_ile
        // Fn_bar = antisymm - fiep[xyz]

        // KM = Re * P_bar' * Ke * P_bar * Re'
        // KGR = -Re * F_bar * G_bar * Re
        // KGP = -Re * G_bar' * Fn_bar' P_bar * Re'
        // Kte = KM + KGR + KGP
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_COROTATION_H_ */
