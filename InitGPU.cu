//---------------------------------------------------------------------------------
//					Getting and then setting the best GPU for the task.
//---------------------------------------------------------------------------------
#include "InitGPU.cuh"
//#define STREAMS 1
using namespace megamol;
using namespace deformables;

inline void check_cuda_errors(const char *filename, const int line_number)
{

	//cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("CUDA error at %s:%i: %s\n", filename, line_number, cudaGetErrorString(error));
		exit(-1);
	}

}


__global__ void computeTheSmallLoop(glm::vec3 *vel,float *m,
									float *roh,glm::vec3 *pos,
									float h,int sphereCount1,float epsilon1,float dt1){//see equaqtion 11	

										int index= blockIdx.x*blockDim.x+threadIdx.x;

										int* nxi;
										int count=0;
										for(int i =0; i <sphereCount1;i++){

											if(glm::distance((glm::vec3)pos[index],(glm::vec3)pos[i])<h){

												if(index!=i){

													nxi[count]=i;
													count++;
												}
											}	

										}
										float *W;
										
										for(int j=0;j<count;j++)
										{
											//Calculating W and gradW

											float q=glm::distance(pos[index],pos[nxi[j]])/h;
											float res=0;
											if(q<=2||q>=0)
											{
												res=0.66+((9/8)*q*q)+((19/24)*q*q*q)-((5/32)*q*q*q*q);
											}

											W[j]=((315/208)*3.14*h*h*h)* res;	
										}
										
										glm::vec3 sum(0.0f,0.0f,0.0f);
										for(int j=0;j<count;j++)
										{
											sum+=(glm::vec3)(((2*(float)m[nxi[j]])/((float)roh[index]+(float)roh[nxi[j]]))*((glm::vec3)vel[nxi[j]]-(glm::vec3)vel[index])*(float)W[j]);
										}		

										vel[index]=vel[index]+(glm::vec3)(epsilon1*sum);
										pos[index]=pos[index]+vel[index]*dt1;

}
__global__ void computeTheBigLoop(float *d_m,float *d_roh,float *d_proh,  float *d_droh,glm::mat3 *d_D,
								  float *d_d,float*d_eta,glm::vec3 *d_pos, glm::vec3 *d_vel,
								  glm::mat3 *d_gradV,glm::vec3 *d_acc, glm::vec3 *d_pacc,glm::mat3 *d_S,
								  glm::vec3 *d_gradS,float *d_P, glm::vec3 *d_gradP,
								  glm::vec3 *d_pie,float h,int sphereCount1,float epsilon1,
								  float c1,float roh01,float alpha1,float jumpN1,float dt1,float etaMax1){

									  int index= blockIdx.x*blockDim.x+threadIdx.x;
									  glm::vec3 *d_gradW;
									  float *d_W;
									  //call other kernels for this one
									  // a kernel will be spawned for each of the material points.

									  // if using dynamic programming then use the following
									  //int threadsPerBlock=32;
									  //int blocksPerGrid=(m.size()/threadsPerBlock) +1;
									  //computeFindNeighbours<<<blocksPerGrid,threadsPerBloack>>>(i,d_pos,d_nxi,d_h);

									  int* nxi=0;
									  int count=0;
									  for(int i =0; i <sphereCount1;i++){

										  if(glm::distance((glm::vec3)d_pos[index],(glm::vec3)d_pos[i])<h){

											  if(index!=i){

												  nxi[count]=i;
												  count++;
											  }
										  }	

									  }
									  //  std::cout<<"thread id "<<index<<" has "<<count<<std::endl;
									  float sumRO=0;
						
									  for(int j=0;j<count;j++)
									  {
										  //Calculating W and gradW

										  float q=glm::distance(d_pos[index],d_pos[nxi[j]])/h;
										  float res=0;
										  if(q<=2||q>=0)
										  {
											  res=0.66+((9/8)*q*q)+((19/24)*q*q*q)-((5/32)*q*q*q*q);
										  }

										  d_W[j]=((315/208)*3.14*h*h*h)* res;	
										  float dXiXj=glm::distance(d_pos[index],d_pos[nxi[j]]);
										  float multiplier=-(9/(4*h))+(19/(8*h*h))*dXiXj+(5/(8*h*h*h))*dXiXj;
										  d_gradW[j].x= multiplier*d_pos[index].x-d_pos[nxi[j]].x;
										  d_gradW[j].y= multiplier*d_pos[index].y-d_pos[nxi[j]].y;
										  d_gradW[j].z= multiplier*d_pos[index].z-d_pos[nxi[j]].z;

										  //calculating gradV
										  glm::vec3 v=(d_vel[nxi[j]]-d_vel[index]);
										  d_gradV[index]=((d_m[nxi[j]]/d_roh[nxi[j]]))*glm::outerProduct(v,d_gradW[j]);

										  d_D[index]=d_gradV[index]+glm::transpose(d_gradV[index]);

										  glm::mat3 D1= 0.5*d_D[index];
										  d_droh[index] = D1[0][0]+D1[1][1]+D1[2][2];

										  //used to calculate density outside the loop

										  sumRO+=d_m[nxi[j]]*d_W[j];
									  }
		
									  //calculating the density roh
									  d_proh[index]=d_roh[index];
									  d_roh[index]=sumRO;

									  //calculating pressure
									  d_P[index]=c1*c1*(d_roh[index]-roh01);

									  glm::mat3 D11= d_D[index];

									  //compute the d from the current rate of deformation tensor
									  float d1=std::sqrt(0.5*(D11[0][0]+D11[1][1]+D11[2][2])*(D11[0][0]+D11[1][1]+D11[2][2]));
									  if(d1==0)
										  d1=1;
									  float exp=std::exp(-(d1*(jumpN1+1)));
									  // since we have fixed n to be .5 other wise use the commented version.
									  float temp=(1-exp)*((1/d1)+(1/d1));
									  //eta[index]=(1-exp)*(std::pow(d1,n-1)*(1/d1));
									  if(temp<etaMax1)
										  d_eta[index]=temp;
									  else
										  d_eta[index]=etaMax1;

									  d_S[index]=d_eta[index]*d_D[index];

									  glm::vec3 sumGP(0.0f,0.0f,0.0f);
									  float ptemp;
						
									  glm::vec3 sumPI(0.0f,0.0f,0.0f);
									  glm::vec3 sumGS(0.0f,0.0f,0.0f);
									  for(int j=0;j<count;j++)
									  {

										  if(glm::dot((d_vel[index]-d_vel[nxi[j]]),(d_pos[index]-d_pos[nxi[j]]))<0)
										  {
											  float dist=glm::distance(d_pos[index],d_pos[nxi[j]]);
											  float mu=h*glm::dot((d_vel[index]-d_vel[nxi[j]]),(d_pos[index]-d_pos[nxi[j]]))/(dist*dist+0.01*h*h);
											  sumPI+=d_m[nxi[j]]*((2*alpha1*c1*(h*mu))/(d_roh[index]+d_roh[j]))*d_gradW[j];
										  }
										  else
											  sumPI+=0;


										  //using the formulation in paper.
										  glm::mat3 sv=(d_m[nxi[j]]/(d_roh[nxi[j]]*d_roh[index]))*(d_S[nxi[j]]+d_S[index]);
										  sumGS+=sv*d_gradW[j];	


										  ptemp=(d_m[nxi[j]])*((d_P[nxi[j]]/std::pow(d_roh[nxi[j]],2))+d_P[index]/std::pow(d_roh[index],2));
										  sumGP+=ptemp*d_gradW[j];

									  }
							
									  d_gradP[index]=sumGP;
									  d_gradS[index]=sumGS;
									  d_pie[index]=sumPI;

									  //updating acceleration
									  glm::vec3 gr(0.0f,-10.0f,0.0f);
									  d_pacc[index]=d_acc[index];
									  d_acc[index]=d_gradS[index]-d_gradP[index]+gr-d_pie[index];

									  //lfs update
									  d_vel[index]=d_vel[index]+0.5f*(d_pacc[index]+d_acc[index])*dt1;
									  // roh has been updated in the paper using LFS which is not reasonable to me 
									  // as the density should be dependent only on the current configuration but it seems 
									  // it needs some sort of help from previous stages.

									  d_roh[index]=d_roh[index]+0.5f*(d_proh[index]+d_roh[index])*dt1;

						
}


void  GPUimplementation::initializeGPU(){

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

void GPUimplementation::theBigLoop(std::vector<float> &m,std::vector<float> &roh,
								   std::vector<float> &proh,std::vector<float> &droh,std::vector<glm::mat3> &D,
								   std::vector<float> &d,std::vector<float>&eta,std::vector<glm::vec3> &pos,
								   std::vector<glm::vec3> &vel, std::vector<glm::mat3> &gradV,std::vector<glm::vec3> &acc,
								   std::vector<glm::vec3> &pacc,std::vector<glm::mat3> &S,std::vector<glm::vec3> &gradS,
								   std::vector<float> &P, std::vector<glm::vec3> &gradP,std::vector<glm::vec3> &gradW,
								   std::vector<glm::vec3> &pie,std::vector<float> &W,float h,int sphereCount,float epsilon,
								   float c,float roh0,float alpha,float jumpN,float dt,float etaMax){

#ifdef STREAMS
									   //streams can be used to copy multiple data structures to the device memeory.
									   cudaStream_t s1;
									   cudaStream_t s2;
									   cudaStream_t s3;

									   cudaStreamCreate(&s1);
									   cudaStreamCreate(&s2);
									   cudaStreamCreate(&s3);
#endif
									   //Declaration for the device storage.
									   float* d_m=0;
									   float* d_roh=0;
									   float* d_proh=0;
									   float* d_droh=0;
									   glm::mat3* d_D=0;
									   float* d_d=0;
									   float*d_eta=0;
									   glm::vec3* d_pos=0;
									   glm::vec3 * d_vel=0;
									   glm::mat3 * d_gradV=0;
									   glm::vec3 * d_acc=0;
									   glm::vec3 * d_pacc=0;
									   glm::mat3 * d_S=0;
									   glm::vec3 * d_gradS=0;
									   float * d_P=0;
									   glm::vec3 * d_gradP=0;
									   glm::vec3 * d_pie=0;

									   //thrust::device_vector <float> d_cm(m.size());
									   //d_cm=m;
									   //float* d_cma=thrust::raw_pointer_cast(&d_cm[0]);


									   //try making a init function where storage for constant
									   //things like m are not transfered again and again.
									   cudaError_t e1 ;
									   cudaError_t e0 =  cudaMalloc((void **)&d_m,sizeof(float)*sphereCount);

#ifndef STREAMS
									   e1 =cudaMemcpy(d_m, &m,sizeof(float)*sphereCount, cudaMemcpyHostToDevice);

#else
									   e1 =cudaMemcpyAsync(d_m, &m,sizeof(m), cudaMemcpyHostToDevice, s1);
#endif
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_roh,sizeof(roh));
#ifndef STREAMS
									   e1=cudaMemcpy(d_roh, &roh, sizeof(roh), cudaMemcpyHostToDevice);
#else
									   e1=cudaMemcpyAsync(d_roh, &roh, sizeof(roh), cudaMemcpyHostToDevice, s2);
#endif
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_proh,sizeof(proh));
#ifndef STREAMS
									   e1=cudaMemcpy(d_proh, &proh, sizeof(proh), cudaMemcpyHostToDevice);
#else
									   e1=cudaMemcpyAsync(d_proh, &proh, sizeof(proh), cudaMemcpyHostToDevice, s3);
#endif
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_droh,sizeof(droh));
									   e1=cudaMemcpy(d_droh, &droh, sizeof(droh), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_D,sizeof(D));
									   e1=cudaMemcpy(d_D, &D, sizeof(D), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;


									   e0=cudaMalloc(&d_d,sizeof(d));
									   e1=cudaMemcpy(d_d, &d, sizeof(d), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_eta,sizeof(eta));
									   e1=cudaMemcpy(d_eta, &eta, sizeof(eta), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_pos,sizeof(pos));
									   e1=cudaMemcpy(d_pos, &pos, sizeof(pos), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_vel,sizeof(vel));
									   e1=cudaMemcpy(d_vel, &vel, sizeof(vel), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_gradV,sizeof(gradV));
									   e1=cudaMemcpy(d_gradV, &gradV, sizeof(gradV), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_acc,sizeof(acc));
									   e1=cudaMemcpy(d_acc, &acc, sizeof(acc), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_pacc,sizeof(pacc));
									   e1=cudaMemcpy(d_pacc, &pacc, sizeof(pacc), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   e0=cudaMalloc(&d_S,sizeof(S));
									   e1=cudaMemcpy(d_S, &S, sizeof(S), cudaMemcpyHostToDevice);
									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$malloc sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
									   if (e1 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   cudaMalloc(&d_gradS,sizeof(gradS));
									   cudaMemcpy(d_gradS, &gradS, sizeof(gradS), cudaMemcpyHostToDevice);

									   cudaMalloc(&d_P,sizeof(P));
									   cudaMemcpy(d_P, &P, sizeof(P), cudaMemcpyHostToDevice);

									   cudaMalloc(&d_gradP,sizeof(gradP));
									   cudaMemcpy(d_gradP, &gradP, sizeof(gradP), cudaMemcpyHostToDevice);

									   cudaMalloc(&d_pie,sizeof(pie));
									   cudaMemcpy(d_pie, &pie, sizeof(pie), cudaMemcpyHostToDevice);




									   int threadsPerBlock=32;
									   int blocksPerGrid=(m.size()/threadsPerBlock) +1;

									   // calling the kernels

									   computeTheBigLoop<<<blocksPerGrid,threadsPerBlock>>>(d_m,d_roh,d_proh,d_droh,
										   d_D, d_d,d_eta,d_pos,d_vel,d_gradV,d_acc,d_pacc,d_S, d_gradS,d_P, d_gradP,
										   d_pie,h, sphereCount,epsilon,c,roh0,alpha,jumpN,dt,etaMax);

									   check_cuda_errors(__FILE__, __LINE__);

									   computeTheSmallLoop<<<blocksPerGrid,threadsPerBlock>>>(d_vel,
										   d_m,d_roh,d_pos, h,sphereCount,epsilon,dt);


									   check_cuda_errors(__FILE__, __LINE__);

									   //copying data back to memory and freeing the device storage.

									   e0=cudaMemcpy(&m,d_m,sizeof(float)*sphereCount,cudaMemcpyDeviceToHost);
									   check_cuda_errors(__FILE__, __LINE__);
									   e1=cudaFree(d_m);

									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;


									   cudaMemcpy( &roh, d_roh, sizeof(d_roh), cudaMemcpyDeviceToHost);
									   cudaFree( d_roh);	

									   cudaMemcpy( &proh,d_proh, sizeof(proh), cudaMemcpyDeviceToHost);
									   cudaFree(d_proh);

									   cudaMemcpy( &droh, d_droh,sizeof(droh), cudaMemcpyDeviceToHost);
									   cudaFree(d_droh);

									   cudaMemcpy( &D,d_D, sizeof(D), cudaMemcpyDeviceToHost);
									   cudaFree(d_D);

									   cudaMemcpy( &d, d_d,sizeof(d), cudaMemcpyDeviceToHost);
									   cudaFree(d_d);

									   cudaMemcpy( &eta, d_eta,sizeof(eta), cudaMemcpyDeviceToHost);
									   cudaFree(d_eta);

									   cudaMemcpy( &pos,d_pos, sizeof(pos), cudaMemcpyDeviceToHost);
									   cudaFree(d_pos);

									   cudaMemcpy( &vel, d_vel,sizeof(vel), cudaMemcpyDeviceToHost);
									   cudaFree(d_vel);

									   cudaMemcpy( &gradV,d_gradV, sizeof(gradV), cudaMemcpyDeviceToHost);
									   cudaFree(d_gradV);

									   cudaMemcpy( &acc,d_acc, sizeof(acc), cudaMemcpyDeviceToHost);
									   cudaFree(d_acc);

									   cudaMemcpy( &pacc, d_pacc,sizeof(pacc), cudaMemcpyDeviceToHost);
									   cudaFree(d_pacc);

									   cudaMemcpy( &S, d_S,sizeof(S), cudaMemcpyDeviceToHost);
									   cudaFree(d_S);

									   e0=cudaMemcpy(&gradS,d_gradS,sizeof(gradS),cudaMemcpyDeviceToHost);
									   e1=cudaFree(d_gradS);

									   if (e0 != cudaSuccess)
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy failed $$$$$$$$$$$$$$$"<< std::endl<< std::endl;

									   else 
										   std::cout<< std::endl<< std::endl<<"$$$$$$$$$$memcpy sucessfully DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;


									   cudaMemcpy(&P,d_P,  sizeof(P), cudaMemcpyDeviceToHost);
									   cudaFree(d_P);

									   cudaMemcpy( &gradP, d_gradP,sizeof(gradP), cudaMemcpyDeviceToHost);
									   cudaFree(d_gradP);

									   cudaMemcpy( &pie, d_pie,sizeof(pie), cudaMemcpyDeviceToHost);
									   cudaFree(d_pie);


}



//__global__ void computeFindNeighbours(int index, std::vector<glm::vec3> &pos,
//									std::vector<glm::vec3> &nxi)){
//
//
//										dim3 i= blockIdx.x*blockDim.x+threadIdx.x;
//
//										if(glm::distance(pos[index],pos[i])<2*h){
//
//											if(index!=i){
//
//												nxi.push_back(i);
//
//											}
//
//
//										}
//								}
//
//
//
//
//								__global__ void calcGradW(W,gradW,pos,i,nxi){
//
//								}
//								__global__ void updategradV(m,roh,vel,gradV,gradW,i,nxi){
//
//
//								}
//
//								__global__ void updateD(gradV,D,i){
//
//
//								}
//
//								__global__ void updateDDE(D,droh)//see equaqtion 10
//								{
//
//
//								}
//
//								__global__ void compute_roh(roh,proh,m,W,i,nxi){
//
//								}
//
//								__global__ void updateP(P,roh0,roh,c)
//									//see equaqtion 6
//								{
//
//
//								}
//
//
//								__global__ void updateEta(eta,D,jn,n)
//									//see equaqtion 4
//								{
//
//								}
//
//								__global__ void updateS(eta,D,S){
//
//
//								}
//
//								//for artificial viscosity equation 12
//								__global__ void compute_pi(m,roh,pie,h,gradW,i,nxi,vel,pos,alpha,c){
//
//								}
//
//								// to update acceleration we need to calculate grad p and grad s
//								__global__ void compute_gradS(roh,S,gradW,i,m,gradS,nxi){
//
//								}
//								__global__ void compute_gradP(m,roh,P,gradW,gradP,i,nxi){
//
//
//								}
