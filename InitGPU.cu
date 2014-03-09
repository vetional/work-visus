//---------------------------------------------------------------------------------
//					Getting and then setting the best GPU for the task.
//---------------------------------------------------------------------------------
#include "InitGPU.cuh"


void   GPUimplementation::initializeGPU(){
    const int kb = 1024;
    const int mb = kb * kb;
	int num_devices, device;
	cudaGetDeviceCount(&num_devices);
	if (num_devices > 1) {
      int max_multiprocessors = 0, max_device = 0;
	   std::cout << "Deformable module for MEGAMOL" << std::endl << "=========" << std::endl << std::endl;

    std::cout << "CUDA version:   v" << CUDART_VERSION << std::endl;    
   std::cout << "Thrust version: v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << std::endl << std::endl; 

	  std::cout<<"CUDA Devices: " << std::endl << std::endl;
      for (device = 0; device < num_devices; device++) {
              cudaDeviceProp props;
              cudaGetDeviceProperties(&props, device);
			  std::cout << device << ": " << props.name << ": " << props.major << "." << props.minor << std::endl;
						std::cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << std::endl;
						std::cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << std::endl;
						std::cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << std::endl;
						std::cout << "  Block registers: " << props.regsPerBlock << std::endl << std::endl;
						std::cout << "  Warp size:         " << props.warpSize << std::endl;
						std::cout << "  Threads per block: " << props.maxThreadsPerBlock << std::endl;
						std::cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]" << std::endl;
						std::cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]" << std::endl;
						std::cout << std::endl;
              if (max_multiprocessors < props.multiProcessorCount) {
                      max_multiprocessors = props.multiProcessorCount;
                      max_device = device;
					  
              }
      }
      cudaSetDevice(max_device);
	  std::cout<<"Selected Device"<<max_device<<std::endl;

	}
}