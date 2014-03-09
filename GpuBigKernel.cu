

#include "GpuBigKernel.h"


void   initializeGPU(){
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
void theBigLoop1(std::vector<float> &m,std::vector<float> &roh,std::vector<float> &proh,
							 std::vector<float> &droh,std::vector<glm::mat3> &D,std::vector<float> &d,
							 std::vector<float>&eta,std::vector<glm::vec3> &pos,
							 std::vector<glm::vec3> &vel,
							 std::vector<glm::mat3> &gradV,std::vector<glm::vec3> &acc,
							 std::vector<glm::vec3> &pacc,std::vector<glm::mat3> &S,
							 std::vector<glm::vec3> &gradS,std::vector<float> &P,
							 std::vector<glm::vec3> &gradP,std::vector<glm::vec3> &gradW,
							 std::vector<glm::vec3> &pie,std::vector<float> &W){


//streams can be used to copy multiple data structures to the device memeory.
		//cudaStream_t s1;
		//cudaStream_t s2;
		//cudaStream_t s3;

		//cudaStreamCreate(&s1);
		//cudaStreamCreate(&s2);
		//cudaStreamCreate(&s3);

//Declare the device storage.
		std::vector<float> d_m;
		std::vector<float> d_roh;
		std::vector<float> d_proh;
		std::vector<float> d_droh;
		std::vector<glm::mat3> d_D;
		std::vector<float> d_d;
		std::vector<float>d_eta;
		std::vector<glm::vec3> d_pos,
		std::vector<glm::vec3> d_vel;
		std::vector<glm::mat3> d_gradV;
		std::vector<glm::vec3> d_acc;
		std::vector<glm::vec3> d_pacc;
		std::vector<glm::mat3> d_S;
		std::vector<glm::vec3> d_gradS;
		std::vector<float> d_P;
		std::vector<glm::vec3> d_gradP;
		std::vector<glm::vec3> d_gradW;
		std::vector<glm::vec3> d_pie;
		std::vector<float> d_W;

				
		//try making a init function where storage for constant 
		//things like m are not transfered again and again.
		
		cudaMalloc(&d_m,sizeof(m));
		cudaMemcpy(d_m, m,sizeof(m), cudaMemcpyHostToDevice);
		//cudaMemcpyAsync(d_elms_cm, elms_cm,size_elms_cm, cudaMemcpyHostToDevice, s1);
		
		cudaMalloc(&d_roh,sizeof(roh));
		cudaMemcpy(d_roh, roh, sizeof(roh), cudaMemcpyHostToDevice);
		//cudaMemcpyAsync(d_vertices_cm, vertices_cm, sizeof(vertices_cm), cudaMemcpyHostToDevice, s2);

		cudaMalloc(&d_proh,sizeof(proh));
		cudaMemcpy(d_proh, proh, sizeof(proh), cudaMemcpyHostToDevice);
		//cudaMemcpyAsync(d_out_elemInfo, out_elemInfo, sizeof(out_elemInfo), cudaMemcpyHostToDevice, s3);

		int threadsPerBlock=128;
		int blocksPerGrid=(m.size()/threadsPerBlock) +1;


		computeTheBigLoop<<<blocksPerGrid,threadsPerBloack>>>();


}

__global__ void computeTheBigLoop(){

	dim3 i= blockIdx.x*blockDim.x+threadIdx.x;
  

  
}
