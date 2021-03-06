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
void theBigLoop1(std::vector<float> &m,std::vector<float> &roh,
		std::vector<float> &proh,std::vector<float> &droh,std::vector<glm::mat3> &D,
		std::vector<float> &d,std::vector<float>&eta,std::vector<glm::vec3> &pos,
		std::vector<glm::vec3> &vel, std::vector<glm::mat3> &gradV,std::vector<glm::vec3> &acc,
		std::vector<glm::vec3> &pacc,std::vector<glm::mat3> &S,std::vector<glm::vec3> &gradS,
		std::vector<float> &P, std::vector<glm::vec3> &gradP,std::vector<glm::vec3> &gradW,
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
		
		
		cudaMalloc(&d_droh,sizeof(droh));
		cudaMemcpy(d_droh, droh, sizeof(droh), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_D,sizeof(D));
		cudaMemcpy(d_D, D, sizeof(D), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_d,sizeof(d));
		cudaMemcpy(d_d, d, sizeof(d), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_eta,sizeof(eta));
		cudaMemcpy(d_eta, eta, sizeof(eta), cudaMemcpyHostToDevice);
		
		cudaMalloc(&d_pos,sizeof(pos));
		cudaMemcpy(d_pos, pos, sizeof(pos), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_vel,sizeof(vel));
		cudaMemcpy(d_vel, vel, sizeof(vel), cudaMemcpyHostToDevice);
		
		
		
		cudaMalloc(&d_gradV,sizeof(gradV));
		cudaMemcpy(d_gradV, gradV, sizeof(gradV), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_acc,sizeof(acc));
		cudaMemcpy(d_acc, acc, sizeof(acc), cudaMemcpyHostToDevice);
		
		cudaMalloc(&d_pacc,sizeof(pacc));
		cudaMemcpy(d_pacc, pacc, sizeof(pacc), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_S,sizeof(S));
		cudaMemcpy(d_S, S, sizeof(S), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_gradS,sizeof(gradS));
		cudaMemcpy(d_gradS, gradS, sizeof(gradS), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_P,sizeof(P));
		cudaMemcpy(d_P, P, sizeof(P), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_gradP,sizeof(gradP));
		cudaMemcpy(d_gradP, gradP, sizeof(gradP), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_gradW,sizeof(gradW));
		cudaMemcpy(d_gradW, gradW, sizeof(gradW), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_pie,sizeof(pie));
		cudaMemcpy(d_pie, pie, sizeof(pie), cudaMemcpyHostToDevice);
		
		
		cudaMalloc(&d_W,sizeof(W));
		cudaMemcpy(d_W, W, sizeof(W), cudaMemcpyHostToDevice);
		


		int threadsPerBlock=128;
		int blocksPerGrid=(m.size()/threadsPerBlock) +1;


		computeTheBigLoop<<<blocksPerGrid,threadsPerBloack>>>(std::vector<float> d_m,d_roh,d_proh,d_droh,
				 d_D, d_d,d_eta,d_pos,d_vel,d_gradV,d_acc,d_pacc,d_S, d_gradS,d_P, d_gradP,
				 d_gradW,d_pie,d_W;);


}

__global__ void computeTheBigLoop(std::vector<float> &m,std::vector<float> &d_roh,std::vector<float> &d_proh,
				 std::vector<float> &d_droh,std::vector<glm::mat3> &d_D,
				 std::vector<float> &d_d,
				 std::vector<float>&d_eta,std::vector<glm::vec3> &d_pos,
				 std::vector<glm::vec3> &d_vel,
				 std::vector<glm::mat3> &d_gradV,std::vector<glm::vec3> &d_acc,
				 std::vector<glm::vec3> &d_pacc,std::vector<glm::mat3> &d_S,
				 std::vector<glm::vec3> &d_gradS,std::vector<float> &d_P,
				 std::vector<glm::vec3> &d_gradP,std::vector<glm::vec3> &d_gradW,
				 std::vector<glm::vec3> &d_pie,std::vector<float> &d_W){

	dim3 index= blockIdx.x*blockDim.x+threadIdx.x;
  
  	//call other kernels for this one
  	// a kernel will be spawned for each of the material points.
  	
  	
  	std::vector<glm::vec3> nxi;
  	
  	
  	// if using dynamic programming then use the following 
	//int threadsPerBlock=32;
	//int blocksPerGrid=(m.size()/threadsPerBlock) +1;
	//computeFindNeighbours<<<blocksPerGrid,threadsPerBloack>>>(i,d_pos,d_nxi,d_h);

	for(int i =0; i <pos.size();i++){
		
		if(glm::distance(pos[index],pos[i])<2*h){

			if(index!=i){
			
				nxi.push_back(i);
			
			}
	
	
		}	

	}
	
	
	for(int j=0;j<nxi.size();j++)
	{
		//Calculating W and  gradW
		d_W[j]=((315/208)*3.14*h*h*h)* w(glm::distance(d_pos[i],d_pos[nxi[j]])/h);	
		float dXiXj=glm::distance(d_pos[i],d_pos[nxi[j]]);
		float multiplier=-(9/(4*h))+(19/(8*h*h))*dXiXj+(5/(8*h*h*h))*dXiXj;
		d_gradW[j].x= multiplier*d_pos[i].x-d_pos[nxi[j]].x;
		d_gradW[j].y= multiplier*d_pos[i].y-d_pos[nxi[j]].y;
		d_gradW[j].z= multiplier*d_pos[i].y-d_pos[nxi[j]].z;
		
		//calculating gradV
		glm::vec3 v=(d_m[nxi[j]]/d_roh[nxi[j]])*(d_vel[nxi[j]]-d_vel[index]);
		d_gradV[index]=glm::outerProduct(v,d_gradW[j]);
	
		d_D[index]=d_gradV[index]+glm::transpose(d_gradV[index]);

		for(int i=0;i<d_D.size();i++){
			glm::mat3 D1= 0.5*d_D[i];
			d_droh[i] = D1[0][0]+D1[1][1]+D1[2][2];

		}
		//used to calculate density outside the loop
		float sum=0;
		sum+=d_m[nxi[j]]*d_W[j];
	}
	
	//calculating the density roh
	d_proh[index]=d_roh[index];
	d_roh[index]=sum;
		
	//calculating pressure 
	d_P[index]=c*c*(d_roh[index]-roh0);
	  
	glm::mat3 D11= 0.5*d_D[index];

	//compute the d from the current rate of deformation tensor 
	float d1=std::sqrt(D11[0][0]+D11[1][1]+D11[2][2]);

	float exp=std::exp(-d1*(jn+1));
	// since we have fixed n to be .5 other wise use the commented version.
	d_eta[index]=(1-exp)*((1/d1)*(1/d1));
	//eta[index]=(1-exp)*(std::pow(d1,n-1)*(1/d1));

	d_S[index]=d_eta[index]*d_D[index];
	
	glm::vec3 sumGP(0.0f,0.0f,0.0f);
	float ptemp;

	glm::vec3 sumPI(0.0f,0.0f,0.0f);
	glm::vec3 sumGS(0.0f,0.0f,0.0f);
	for(int j=0;j<nxi.size();j++)
	{

		if(glm::dot((d_vel[index]-d_vel[nxi[j]]),(d_pos[index]-d_pos[nxi[j]]))<0)
		{
			float dist=glm::distance(d_pos[index],d_pos[nxi[j]]);
			float mu=h*glm::dot((d_vel[index]-d_vel[nxi[j]]),(d_pos[index]-d_pos[nxi[j]]))/(dist*dist+0.01*h*h);
			sumPI+=d_m[nxi[j]]*((2*alpha*c*(h*mu))/(d_roh[index]+d_roh[j]))*d_gradW[j];
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



}

void SPHSimulation::updateS(std::vector<float> &eta,std::vector<glm::mat3>D,
	std::vector<glm::mat3>&S){


	for(int i=0;i<D.size();i++){

		S[i]=eta[i]*D[i];


}








__global__ void computeFindNeighbours(int index, std::vector<glm::vec3> &pos, 
						std::vector<glm::vec3> &nxi)){

	
	dim3 i= blockIdx.x*blockDim.x+threadIdx.x;

	if(glm::distance(pos[index],pos[i])<2*h){

		if(index!=i){
		
			nxi.push_back(i);
			
		}
	
	
	}
}
	



__global__ void calcGradW(W,gradW,pos,i,nxi){

}
__global__ void updategradV(m,roh,vel,gradV,gradW,i,nxi){


}

__global__ void updateD(gradV,D,i){


}

__global__ void updateDDE(D,droh)//see equaqtion 10
{


}

__global__ void compute_roh(roh,proh,m,W,i,nxi){

}

__global__ void updateP(P,roh0,roh,c)
//see equaqtion 6
{


}


__global__ void updateEta(eta,D,jn,n)
//see equaqtion 4
{

}

__global__ void updateS(eta,D,S){


}

//for artificial viscosity equation 12
__global__ void compute_pi(m,roh,pie,h,gradW,i,nxi,vel,pos,alpha,c){

}

// to update acceleration we need to calculate grad p and grad s
__global__ void compute_gradS(roh,S,gradW,i,m,gradS,nxi){

}
__global__ void compute_gradP(m,roh,P,gradW,gradP,i,nxi){


}
