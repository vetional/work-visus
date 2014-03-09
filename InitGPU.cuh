#ifndef 	INITGPU
#define INITGPU
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/version.h>


 //void initializeGPU();

namespace megamol {
	namespace deformables {
		class GPUimplementation
		{
			public:
				 void  GPUimplementation::initializeGPU();

		};
	} /* end namespace deformables */
} /* end namespace megamol */
#endif