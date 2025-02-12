#ifndef SRC_MORPHING_CYCLIC_GRAPH_DECOMPOSITION_H_
#define SRC_MORPHING_CYCLIC_GRAPH_DECOMPOSITION_H_

#include <vector>

namespace espreso{
    

class CyclicGraphDecomposition{
//TODO implement a function to calculate CDCs for unknown values
public:
    
    /*
        @block_distribution[c + r*nranks] denotes the index of an MPI process assembling that block
        
        if @active_block_indices[i] is true, then this MPI process requires data from the i-th MPI process
        
        @block_processes[i] denotes a list of MPI processes requiring data from the i-th MPI process
        
        @process_blocks[i] denotes a list of blocks required by the i-th MPI process
    */
    
    static void define_workload(
        std::vector<esint> &block_distribution,
        std::vector<bool> &active_block_indices,
        std::vector<std::vector<esint>> &block_processes,
        std::vector<std::vector<esint>> &process_blocks,
        esint nranks,
        esint proc_idx
    );
    
    
};
    

} /* end of namespace espreso */

#endif /* SRC_MORPHING_CYCLIC_GRAPH_DECOMPOSITION_H_ */
