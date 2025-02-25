
#include "allocator_new.h"



namespace espreso {



AllocatorDummy_new AllocatorDummy_new::singleton_ff = AllocatorDummy_new(false, false);
AllocatorDummy_new AllocatorDummy_new::singleton_ft = AllocatorDummy_new(false, true);
AllocatorDummy_new AllocatorDummy_new::singleton_tf = AllocatorDummy_new(true,  false);
AllocatorDummy_new AllocatorDummy_new::singleton_tt = AllocatorDummy_new(true,  true);

AllocatorCPU_new AllocatorCPU_new::singleton = AllocatorCPU_new();

AllocatorGPU_new AllocatorGPU_new::singleton = AllocatorGPU_new();

AllocatorHostPinned_new AllocatorHostPinned_new::singleton = AllocatorHostPinned_new();



}
