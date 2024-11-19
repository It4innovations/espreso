
#include "op.print.eigenvalues.h"

#include "math/wrappers/math.lapack.h"

namespace espreso {

void PrintEigenValuesKernel::simd(const char * name, size_t size, SIMD M[], size_t n)
{
    m.resize(size, size);
    e.resize(size);
    for (size_t s = 0; s < SIMD::size; ++s) {
        if (!std::isnan(M[0][s])) {
            for (size_t r = 0; r < size; ++r) {
                for (size_t c = 0; c < size; ++c) {
                    m.vals[r * size + c] = M[r * size + c][s];
                }
            }
            math::lapack::get_eig_sym(m, e);

            printf("%s:", name);
            for (size_t v = 0; v < n; ++v) {
                printf(" %+.5e", e.vals[v]);
            }
            if (n < (size_t)e.size) {
                printf(" ... %+.5e", e.vals[e.size - 1]);
            }
            printf("\n");
        }
    }
}

}
