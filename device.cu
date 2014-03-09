#include <iostream>
#include "thrust/sort.h"
#include "thrust/device_vector.h"
#include "thrust/copy.h"
#include "thrust/version.h"
#include "device.h"
using namespace std;
void sort_on_device(thrust::host_vector<int>& h_vec)
{
//---------------------------------------------------------------------------------
//					Getting and then setting the best GPU for the task.
//---------------------------------------------------------------------------------

    const int kb = 1024;
    const int mb = kb * kb;
	int num_devices, device;
	cudaGetDeviceCount(&num_devices);
	if (num_devices > 1) {
      int max_multiprocessors = 0, max_device = 0;
	   cout << "NBody.GPU" << endl << "=========" << endl << endl;

    cout << "CUDA version:   v" << CUDART_VERSION << endl;    
   cout << "Thrust version: v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << endl << endl; 

	  cout<<"CUDA Devices: " << endl << endl;
      for (device = 0; device < num_devices; device++) {
              cudaDeviceProp props;
              cudaGetDeviceProperties(&props, device);
			  cout << device << ": " << props.name << ": " << props.major << "." << props.minor << endl;
						cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << endl;
						cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << endl;
						cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << endl;
						cout << "  Block registers: " << props.regsPerBlock << endl << endl;

						cout << "  Warp size:         " << props.warpSize << endl;
						cout << "  Threads per block: " << props.maxThreadsPerBlock << endl;
						cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]" << endl;
						cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]" << endl;
						cout << endl;
              if (max_multiprocessors < props.multiProcessorCount) {
                      max_multiprocessors = props.multiProcessorCount;
                      max_device = device;
					  
              }
      }
      cudaSetDevice(max_device);
	  cout<<"Selected Device"<<max_device<<endl;
}
//---------------------------------------------------------------------------------
//				Getting the job done.
//---------------------------------------------------------------------------------
    
	
	// transfer data to the device
		thrust::device_vector<int> d_vec = h_vec;

    // sort data on the device
		 thrust::sort(d_vec.begin(), d_vec.end());
    
    // transfer data back to host
	     thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
}
