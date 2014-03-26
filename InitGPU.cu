//---------------------------------------------------------------------------------
//					Getting and then setting the best GPU for the task.
//---------------------------------------------------------------------------------
#include "InitGPU.cuh"
#include "math.h"
//#define STREAMS 1
using namespace megamol;
using namespace deformables;

float* d_m;
float* d_roh;
float* d_proh;
float* d_droh;
glm::mat3* d_D;
float* d_d;
float*d_eta;
glm::vec3* d_pos;
glm::vec3 * d_vel;
glm::mat3 * d_gradV;
glm::vec3 * d_acc;
glm::vec3 * d_pacc;
glm::mat3 * d_S;
glm::vec3 * d_gradS;
float * d_P;
glm::vec3 * d_gradP;
glm::vec3 * d_pie;

inline void check_cuda_errors(const char *filename, const int line_number,const char *msg)
{

	cudaThreadSynchronize();
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("CUDA error at %s:%i:during  %s >> %s\n", filename, line_number,msg, cudaGetErrorString(error));
		exit(-1);
	}

}


//__device__ void computeTheSmallLoop(glm::vec3 *vel,float *m,
//									float *roh,glm::vec3 *pos,int index,
//									float h,int sphereCount1,float epsilon1,float dt1){//see equaqtion 11	
//
//
//										int nxi[200];
//										int count=0;
//										for(int i =0; i <sphereCount1;i++){
//
//											if(glm::distance((glm::vec3)pos[index],(glm::vec3)pos[i])<h){
//
//												if(index!=i){
//
//													nxi[count]=i;
//													count++;
//												}
//											}	
//
//										}
//										float W[200];
//
//										for(int j=0;j<count;j++)
//										{
//											//Calculating W and gradW
//
//											float q=glm::distance(pos[index],pos[nxi[j]])/h;
//											float res=0;
//											if(q<=2||q>=0)
//											{
//												res=0.66+((9/8)*q*q)+((19/24)*q*q*q)-((5/32)*q*q*q*q);
//											}
//
//											W[j]=((315/208)*3.14*h*h*h)* res;	
//										}
//
//										glm::vec3 sum(0.0f,0.0f,0.0f);
//										for(int j=0;j<count;j++)
//										{
//											float cst=(((2.0*m[nxi[j]])/(roh[index]+roh[nxi[j]]))*W[j]);
//											sum+=(vel[nxi[j]]-vel[index])*cst;
//										}		
//
//										vel[index]=vel[index]+(glm::vec3)(epsilon1*sum);
//										pos[index]=pos[index]+vel[index]*dt1;
//										//pos[index]=pos[index]+glm::vec3(0,0,-0.001);
//
//
//}
//
//__device__ void tbl1(glm::vec3 *d_pos, glm::vec3 *d_vel,glm::mat3 *d_gradV,glm::mat3 *d_D,float *d_m,float *d_roh,float *d_droh,float h,int *nxi,int count,int index){
//
//
//	float sumRO=0;
//	glm::vec3 d_gradW[200];
//	float d_W[200];
//	for(int j=0;j<count;j++)
//	{
//		//Calculating W and gradW
//		int ni=nxi[j];
//		float q1=glm::distance(d_pos[index],d_pos[ni]);
//		float q=q1/h;
//		float res=0;
//		if(q<=2||q>=0)
//		{
//			res=0.66+((9/8)*q*q)+((19/24)*q*q*q)-((5/32)*q*q*q*q);
//		}
//
//		d_W[j]=((315/208)*3.14*h*h*h)* res;	
//
//
//		float multiplier=-(9/(4*h))+(19/(8*h*h))*q1+(5/(8*h*h*h))*q1;
//
//		glm::vec3 temp;
//
//		temp.x= multiplier*d_pos[index].x-d_pos[nxi[j]].x;
//		temp.y= multiplier*d_pos[index].y-d_pos[nxi[j]].y;
//		temp.z= multiplier*d_pos[index].z-d_pos[nxi[j]].z;
//
//		d_gradW[j]=(glm::vec3)temp;
//
//		//calculating gradV
//		glm::vec3 v=(glm::vec3)(d_vel[nxi[j]]-d_vel[index]);
//		d_gradV[index]=(glm::mat3)((d_m[nxi[j]]/d_roh[nxi[j]]))*glm::outerProduct(v,d_gradW[j]);
//
//		d_D[index]=(glm::mat3)(d_gradV[index]+glm::transpose(d_gradV[index]));
//
//		glm::mat3 D1= 0.5f*d_D[index];
//		d_droh[index] = (float)(D1[0][0]+D1[1][1]+D1[2][2]);
//
//		//used to calculate density outside the loop
//
//		sumRO+=d_m[nxi[j]]*d_W[j];
//
//	}
//
//
//}
//__device__ void tbl2(float *d_roh,glm::mat3 *d_D, float *d_P,float*d_eta,glm::mat3 *d_S,int index,float roh01,float c1,float jumpN1,float etaMax1){
//
//	//calculating pressure
//	d_P[index]=c1*c1*(d_roh[index]-roh01);
//
//	glm::mat3 D11= d_D[index];
//
//	//compute the d from the current rate of deformation tensor
//	float d1=__fsqrt_ru(0.5*(D11[0][0]+D11[1][1]+D11[2][2])*(D11[0][0]+D11[1][1]+D11[2][2]));
//	if(d1==0)
//		d1=1;
//	float ep=__expf(-(d1*(jumpN1+1)));
//	// since we have fixed n to be .5 other wise use the commented version.
//	float temp=(1-ep)*((1/d1)+(1/d1));
//	//eta[index]=(1-exp)*(std::pow(d1,n-1)*(1/d1));
//	if(temp<etaMax1)
//		d_eta[index]=temp;
//	else
//	{ 
//		d_eta[index]=temp;
//		etaMax1=temp;
//	}
//	d_S[index]=d_eta[index]*d_D[index];
//
//
//}
//
//__device__ void tbl3(glm::vec3 * d_vel,glm::vec3 *d_pos, float* d_m,glm::mat3 * d_S,float* d_P,glm::vec3 * d_acc,glm::vec3 * d_pacc,
//					 glm::vec3 * d_gradS,glm::vec3 * d_gradP,float* d_roh,float* d_proh,glm::vec3 * d_pie,
//					 int count,int index,int *nxi,float alpha1,float c1,float h,float dt1){
//
//						 glm::vec3 sumGP(0.0f,0.0f,0.0f);
//						 float ptemp;
//						 glm::vec3 d_gradW[200];
//						 float d_W[200];
//						 glm::vec3 sumPI(0.0f,0.0f,0.0f);
//						 glm::vec3 sumGS(0.0f,0.0f,0.0f);
//						 for(int j=0;j<count;j++)
//						 {
//							 int ni=nxi[j];
//							 if(glm::dot((d_vel[index]-d_vel[nxi[j]]),(d_pos[index]-d_pos[nxi[j]]))<0)
//							 {
//								 float dist=glm::distance(d_pos[index],d_pos[nxi[j]]);
//								 float mu=h*glm::dot((d_vel[index]-d_vel[ni]),(d_pos[index]-d_pos[ni]))/(dist*dist+0.01*h*h);
//								 sumPI+=d_m[nxi[j]]*((2*alpha1*c1*(h*mu))/(d_roh[index]+d_roh[j]))*d_gradW[j];
//							 }
//							 else
//								 sumPI+=0;
//
//
//							 //using the formulation in paper.
//							 glm::mat3 sv=(d_m[ni]/(d_roh[ni]*d_roh[index]))*(d_S[ni]+d_S[index]);
//							 sumGS+=sv*d_gradW[j];	
//
//
//							 ptemp=(d_m[ni])*((d_P[ni]/(d_roh[ni]*d_roh[ni]))+d_P[index]/(d_roh[index]*d_roh[index]));
//							 sumGP+=ptemp*d_gradW[j];
//
//						 }
//
//						 d_gradP[index]=sumGP;
//						 d_gradS[index]=sumGS;
//						 d_pie[index]=sumPI;
//
//						 //updating acceleration
//						 glm::vec3 gr(0.0f,-10.0f,0.0f);
//						 d_pacc[index]=d_acc[index];
//						 d_acc[index]=d_gradS[index]-d_gradP[index]+gr-d_pie[index];
//
//						 //lfs update
//						 d_vel[index]=d_vel[index]+0.5f*(d_pacc[index]+d_acc[index])*dt1;
//						 // roh has been updated in the paper using LFS which is not reasonable to me 
//						 // as the density should be dependent only on the current configuration but it seems 
//						 // it needs some sort of help from previous stages.
//
//						 d_roh[index]=d_roh[index]+0.5f*(d_proh[index]+d_roh[index])*dt1;
//}

__global__ void tbl(glm::vec3 *d_pos, glm::vec3 *d_vel,glm::mat3 *d_gradV,glm::mat3 *d_D, float *d_m,float *d_roh,
					float *d_droh,float *d_P,float*d_eta,glm::mat3 *d_S,glm::vec3 * d_acc,glm::vec3 * d_pacc,
					glm::vec3 * d_gradS,glm::vec3 * d_gradP,float* d_proh,glm::vec3 * d_pie,
					float h,int sphereCount,float roh01,float c1,float jumpN1,float etaMax1,float alpha1,float dt1,float epsilon1){

						int index= blockIdx.x*blockDim.x+threadIdx.x;
						if(index<sphereCount){

							int nxi[200];
							int count=0;
							for(int i =0; i <sphereCount;i++)
							{

								float dis=glm::distance(d_pos[index],d_pos[i]);

								if(dis<h)
								{

									if(index!=i)
									{

										nxi[count]=i;
										count++;

									}
								}	

								//d_pos[i].z=5.0f;
							}


							//tbl1(d_pos,d_vel,d_gradV,d_D,d_m,d_roh,d_droh,h,nxi,count,index);
							//tbl2(d_roh,d_D,d_P,d_eta,d_S,index,roh01,c1,jumpN1,etaMax1);
							//tbl3( d_vel,d_pos,d_m,d_S, d_P, d_acc, d_pacc, d_gradS, d_gradP,d_roh,
							//	d_proh,d_pie,count, index, nxi,alpha1, c1, h,dt1);
							//computeTheSmallLoop(d_vel,d_m,d_roh,d_pos,index,h,sphereCount,epsilon1,dt1);
						}
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
void GPUimplementation::MallocCudaptr(float sphereCount){
	check_cuda_errors(__FILE__, __LINE__,"malloc started");
	cudaMalloc((void**)&d_m,sizeof(float)*sphereCount);
	cudaMalloc((void**)&d_roh,sizeof(float)*sphereCount);
	cudaMalloc((void**)&d_proh,sizeof(float)*sphereCount);
	cudaMalloc((void**)&d_droh,sizeof(float)*sphereCount);
	cudaMalloc((void**)&d_D,sizeof(glm::mat3)*sphereCount);
	cudaMalloc((void**)&d_d,sizeof(float)*sphereCount);
	cudaMalloc((void**)&d_eta,sizeof(float)*sphereCount);
	cudaMalloc((void**)&d_pos,sizeof(glm::vec3)*sphereCount);
	cudaMalloc((void**)&d_vel,sizeof(glm::vec3)*sphereCount);
	cudaMalloc((void**)&d_gradV,sizeof(glm::mat3)*sphereCount);
	cudaMalloc((void**)&d_acc,sizeof(glm::vec3)*sphereCount);
	cudaMalloc((void**)&d_pacc,sizeof(glm::vec3)*sphereCount);
	cudaMalloc((void**)&d_S,sizeof(glm::mat3)*sphereCount);
	cudaMalloc((void**)&d_gradS,sizeof(glm::vec3)*sphereCount);
	cudaMalloc((void**)&d_P,sizeof(float)*sphereCount);
	cudaMalloc((void**)&d_gradP,sizeof(glm::vec3)*sphereCount);
	cudaMalloc((void**)&d_pie,sizeof(glm::vec3)*sphereCount);
	check_cuda_errors(__FILE__, __LINE__,"malloc Done");



};
#ifdef eg
void GPUimplementation::theBigLoop(float *m,float *roh,float *proh,float *droh,glm::mat3 *D,
								   float *d,float*eta,glm::vec3 *pos, glm::vec3 *vel, glm::mat3 *gradV,glm::vec3 *acc,
								   glm::vec3 *pacc,glm::mat3 *S,glm::vec3 *gradS,
								   float *P, glm::vec3 *gradP,glm::vec3 *pie,float h,int sphereCount,float epsilon,
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
									   //try making a init function where storage for constant
									   //things like m are not transfered again and again.

#ifndef STREAMS
									   cudaMemcpy(d_m, m,sizeof(float)*sphereCount, cudaMemcpyHostToDevice);

#else
									   cudaMemcpyAsync(d_m, m,sizeof(float)*sphereCount, cudaMemcpyHostToDevice, s1);
#endif
#ifndef STREAMS
									   cudaMemcpy(d_roh, roh, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
#else
									   cudaMemcpyAsync(d_roh, roh,sizeof(float)*sphereCount, cudaMemcpyHostToDevice, s2);
#endif
#ifndef STREAMS
									   cudaMemcpy(d_proh, proh, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
#else
									   cudaMemcpyAsync(d_proh, proh,sizeof(float)*sphereCount, cudaMemcpyHostToDevice, s3);
#endif
									   cudaMemcpy(d_droh, droh, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_D,D, sizeof(glm::mat3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_d, d,sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_eta, eta, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_pos, pos, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_vel, vel, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_gradV, gradV, sizeof(glm::mat3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_acc,acc, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_pacc, pacc, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_S, S, sizeof(glm::mat3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_gradS, gradS,sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_P, P,sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_gradP, gradP,sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_pie, pie, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   check_cuda_errors(__FILE__, __LINE__,"Memcpy from host to device done");
									   int threadsPerBlock=32;
									   int blocksPerGrid=(sphereCount/threadsPerBlock) +1;

									   // calling the kernel

									   tbl<<<blocksPerGrid,threadsPerBlock>>>(d_pos, d_vel,d_gradV,d_D,d_m,d_roh,
										   d_droh,d_P,d_eta,d_S,d_acc,d_pacc, d_gradS,d_gradP,d_proh,d_pie,
										   h,sphereCount, roh0,c,jumpN,etaMax,alpha,dt,epsilon);

									   check_cuda_errors(__FILE__, __LINE__,"After Kernel");

									   //copying data back to memory and freeing the device storage.


									   std::cout<<"value of pos[110].z before transfer from cpu to gpu "<<pos[110].z<<std::endl;

									   cudaMemcpy(&m,d_m,sizeof(float)*sphereCount,cudaMemcpyDeviceToHost);
									   check_cuda_errors(__FILE__, __LINE__,"malloc done m");
									   cudaMemcpy( &roh, d_roh, sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);

									   check_cuda_errors(__FILE__, __LINE__,"malloc done roh ");
									   cudaMemcpy( &proh,d_proh, sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &droh, d_droh,sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);

									   cudaMemcpy( &d, d_d,sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &eta, d_eta,sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &pos,d_pos, sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &vel, d_vel,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &gradV,d_gradV, sizeof(glm::mat3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &acc,d_acc, sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &pacc, d_pacc,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &S, d_S,sizeof(glm::mat3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy(&gradS,d_gradS,sizeof(glm::vec3)*sphereCount,cudaMemcpyDeviceToHost);
									   cudaMemcpy(&P,d_P, sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &gradP, d_gradP,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &pie, d_pie,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   std::cout<<"value of pos[110].z after copy "<<pos[110].z<<std::endl;
									   check_cuda_errors(__FILE__, __LINE__,"malloc done pie");

}
#else

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
									   //try making a init function where storage for constant
									   //things like m are not transfered again and again.

#ifndef STREAMS
									   cudaMemcpy(d_m, &m,sizeof(float)*sphereCount, cudaMemcpyHostToDevice);

#else
									   cudaMemcpyAsync(d_m, &m,sizeof(float)*sphereCount, cudaMemcpyHostToDevice, s1);
#endif
#ifndef STREAMS
									   cudaMemcpy(d_roh, &roh, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
#else
									   cudaMemcpyAsync(d_roh, &roh,sizeof(float)*sphereCount, cudaMemcpyHostToDevice, s2);
#endif
#ifndef STREAMS
									   cudaMemcpy(d_proh, &proh, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
#else
									   cudaMemcpyAsync(d_proh, &proh,sizeof(float)*sphereCount, cudaMemcpyHostToDevice, s3);
#endif
									   cudaMemcpy(d_droh, &droh, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_D, &D, sizeof(glm::mat3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_d, &d,sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_eta, &eta, sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_pos, &pos, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_vel, &vel, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_gradV, &gradV, sizeof(glm::mat3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_acc, &acc, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_pacc, &pacc, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_S, &S, sizeof(glm::mat3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_gradS, &gradS,sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_P, &P,sizeof(float)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_gradP, &gradP,sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   cudaMemcpy(d_pie, &pie, sizeof(glm::vec3)*sphereCount, cudaMemcpyHostToDevice);
									   check_cuda_errors(__FILE__, __LINE__,"Memcpy sucess");
									   int threadsPerBlock=32;
									   int blocksPerGrid=(sphereCount/threadsPerBlock) +1;

									   // calling the kernel

									   tbl<<<blocksPerGrid,threadsPerBlock>>>(d_pos, d_vel,d_gradV,d_D,d_m,d_roh,
										   d_droh,d_P,d_eta,d_S,d_acc,d_pacc, d_gradS,d_gradP,d_proh,d_pie,
										   h,sphereCount, roh0,c,jumpN,etaMax,alpha,dt,epsilon);

									   check_cuda_errors(__FILE__, __LINE__,"kernenl execution sucess");

									   //copying data back to memory and freeing the device storage.


									   std::cout<<"value of pos[110].z before transfer from cpu to gpu "<<pos[110].z<<std::endl;

									   cudaMemcpy(&m,d_m,sizeof(float)*sphereCount,cudaMemcpyDeviceToHost);
									   check_cuda_errors(__FILE__, __LINE__,"m copied back sucess");
									   cudaMemcpy( &roh, d_roh, sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);

									   check_cuda_errors(__FILE__, __LINE__,"roh copied back sucess");
									   cudaMemcpy( &proh,d_proh, sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   droh.resize(sphereCount);
									   cudaMemcpy( &droh, d_droh,sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);

									   cudaMemcpy( &d, d_d,sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &eta, d_eta,sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &pos,d_pos, sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &vel, d_vel,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &gradV,d_gradV, sizeof(glm::mat3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &acc,d_acc, sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &pacc, d_pacc,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &S, d_S,sizeof(glm::mat3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy(&gradS,d_gradS,sizeof(glm::vec3)*sphereCount,cudaMemcpyDeviceToHost);
									   cudaMemcpy(&P,d_P, sizeof(float)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &gradP, d_gradP,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   cudaMemcpy( &pie, d_pie,sizeof(glm::vec3)*sphereCount, cudaMemcpyDeviceToHost);
									   std::cout<<"value of pos[110].z after copy "<<pos[110].z<<std::endl;
									   check_cuda_errors(__FILE__, __LINE__,"all mem copied back sucess");

}

#endif
void GPUimplementation::freeMem(){
	check_cuda_errors(__FILE__, __LINE__,"cuda free started");
	cudaFree(d_m);
	cudaFree( d_roh);	
	cudaFree(d_proh);
	cudaFree(d_droh);
	cudaFree(d_D);
	cudaFree(d_d);
	cudaFree(d_eta);
	cudaFree(d_pos);
	cudaFree(d_vel);
	cudaFree(d_gradV);
	cudaFree(d_acc);
	cudaFree(d_pacc);
	cudaFree(d_S);
	cudaFree(d_gradS);
	cudaFree(d_P);
	cudaFree(d_gradP);
	cudaFree(d_pie);
	check_cuda_errors(__FILE__, __LINE__,"cuda free done" );
}
