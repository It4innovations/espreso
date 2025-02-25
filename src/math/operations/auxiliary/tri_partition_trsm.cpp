
#include "math/operations/auxiliary/tri_partition_trsm.h"



namespace espreso {
namespace math {
namespace operations {



void tri_partition_trsm::init(char algorithm_, char direction_, int parameter_)
{
    algorithm = algorithm_;
    direction = direction_;
    parameter = parameter_;

    set_config_called = true;
}



void tri_partition_trsm::set_system(size_t sys_size_, size_t sys_nrhs_)
{
    sys_size = sys_size_;
    sys_nrhs = sys_nrhs_;

    set_system_called = true;
}



void tri_partition_trsm::setup()
{
    if(!set_config_called) eslog::error("config is not set\n");
    if(!set_system_called) eslog::error("system is not set\n");

    if(direction == 'H') partition_range = sys_nrhs;
    if(direction == 'V') partition_range = sys_size;

    switch(algorithm) {
        case 'U':
        case 'M':
        {
            if(parameter > 0) {
                num_chunks = parameter;
                break;
            }
            if(parameter < 0) {
                size_t approx_chunk_size = -parameter;
                num_chunks = (partition_range - 1) / approx_chunk_size + 1;
                break;
            }
            eslog::error("invalid partition parameter\n");
        }
        default:
        {
            eslog::error("invalid partition algorithm\n");
        }
    }

    setup_called = true;
}



size_t tri_partition_trsm::get_num_chunks()
{
    if(!setup_called) eslog::error("setup was not called\n");

    return num_chunks;
}



void tri_partition_trsm::set_output_partition(VectorDenseView_new<size_t> * partition_)
{
    partition = partition_;
}



void tri_partition_trsm::perform()
{
    if(!setup_called) eslog::error("setup was not called\n");
    if(partition == nullptr) eslog::error("partition is not set\n");
    if(partition->size != num_chunks + 1) eslog::error("wrong partition size\n");
    if(sizeof(size_t) != sizeof(double)) eslog::error("incompatible types\n");

    double * partition_sizet = partition->vals;
    double * partition_double = reinterpret_cast<double*>(partition->vals);

    switch(algorithm) {
        case 'U': // uniform
        {
            for(size_t i = 0; i <= num_chunks; i++) {
                partition_double[i] = (double)i / num_chunks;
            }
            break;
        }
        case 'M': // Minimum possible amount of total work
        {
            // W = sum W[i]
            // W[i] = P[i]^2 * (P[i] - P[i-1]) // P[N]=0, P[0]=n, then need to flip and reverse
            // P = arg min W   <===>   P solves dW/dP[j] = 0 for all j
            // dW/dP[j] = -P[j+1]^2 + 2P[j]*(P[j]-P[j-1]) + P[j]^2   =0
            // P[j+1] = sqrt(P[j] * (3P[j]+P[j-1]))
            // not very rigorous, but it works
            double * p = partition_double;
            p[0] = 0.0;
            p[1] = 1.0;
            for(size_t i = 2; i <= num_chunks; i++) {
                p[i] = std::sqrt(p[i-1] * (3 * p[i-1] - 2 * p[i-2]));
            }
            for(size_t i = 0; i <= num_chunks; i++) {
                p[i] = p.back() - p[i];
            }
            std::reverse(p, p + partition->size);
            for(size_t i = 0; i <= num_chunks; i++) {
                p[i] = p.back() - p[i];
            }
            break;
        }
        default:
        {
            eslog::error("invalid partition algorithm\n");
        }
    }

    double max_bound = partition_double[num_chunks];
    for(size_t i = 0; i <= num_chunks; i++) {
        partition_sizet[i] = (size_t)(partition_range * (partition_double[i] / max_bound));
    }
}



}
}
}
