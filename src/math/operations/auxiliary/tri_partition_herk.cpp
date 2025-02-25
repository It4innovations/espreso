
#include "math/operations/auxiliary/tri_partition_herk.h"



namespace espreso {
namespace math {
namespace operations {



void tri_partition_herk::set_config(char algorithm_, char direction_, int parameter_, char herk_strategy_)
{
    algorithm = algorithm_;
    direction = direction_;
    parameter = parameter_;
    herk_strategy = herk_strategy_;

    set_config_called = true;
}



void tri_partition_herk::set_system(size_t size_n_, size_t size_k_)
{
    size_n = size_n_;
    size_k = size_k_;

    set_system_called = true;
}



void tri_partition_herk::setup()
{
    if(!set_config_called) eslog::error("config is not set\n");
    if(!set_system_called) eslog::error("system is not set\n");

    if(direction == 'N') partition_range = size_n;
    if(direction == 'K') partition_range = size_k;

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
                size_t partition_range = get_partition_range();
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



size_t tri_partition_herk::get_num_chunks()
{
    if(!setup_called) eslog::error("setup was not called\n");

    return num_chunks;
}



void tri_partition_herk::set_output_partition(VectorDenseView_new<size_t> * partition_)
{
    partition = partition_;
}



void tri_partition_herk::perform()
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
            if(herk_strategy == 'T') {
                // in the end, it is the same as uniform
                for(size_t i = 0; i <= num_chunks; i++) {
                    partition_double[i] = (double)i / num_chunks;
                }
                break;
            }
            if(herk_strategy == 'Q') {
                // V = sum V[i]
                // V[i] = P[i]^2 * (P[i] - P[i-1]) // P[0]=0, P[N]=n now
                // P = arg min V   <===>   P solves dV/dP[j] = 0 for all j
                // dV/dP[j] = -P[j+1]^2 + 2P[j]*(P[j]-P[j-1]) + P[j]^2   =0
                // P[j+1] = sqrt(P[j] * (3P[j]-2P[j-1]))
                // not very rigorous, but it works
                double * p = partition_double;
                p[0] = 0.0;
                p[1] = 1.0;
                for(size_t i = 2; i <= num_chunks; i++)
                {
                    p[i] = std::sqrt(p[i-1] * (3 * p[i-1] - 2 * p[i-2]));
                }
                break;
            }
            eslog::error("wrong herk_strategy\n");
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
