#pragma once
#ifndef 	INITGPU
#define INITGPU
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/version.h>
#include <thrust/host_vector.h>

#include "glm-0.9.4.6\glm\glm\glm.hpp"

//#define eg 1


namespace megamol {
	namespace deformables {
		class GPUimplementation
		{
		public:
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

			void GPUimplementation::freeMem();
			void GPUimplementation::MallocCudaptr(float sphereCount);
			float  GPUimplementation::w(float q);
			void  GPUimplementation::initializeGPU();

#ifdef eg
			void GPUimplementation::theBigLoop(float *m,float *roh,float *proh,float *droh,glm::mat3 *D,
				float *d,float*eta,glm::vec3 *pos, glm::vec3 *vel, glm::mat3 *gradV,glm::vec3 *acc,
				glm::vec3 *pacc,glm::mat3 *S,glm::vec3 *gradS,
				float *P, glm::vec3 *gradP,glm::vec3 *pie,float h,int sphereCount,float epsilon,
				float c,float roh0,float alpha,float jumpN,float dt,float etaMax);
#else
			void  GPUimplementation::theBigLoop(std::vector<float> &m,std::vector<float> &roh,std::vector<float> &proh,
				std::vector<float> &droh,std::vector<glm::mat3> &D,std::vector<float> &d,
				std::vector<float>&eta,std::vector<glm::vec3> &pos,
				std::vector<glm::vec3> &vel,
				std::vector<glm::mat3> &gradV,std::vector<glm::vec3> &acc,
				std::vector<glm::vec3> &pacc,std::vector<glm::mat3> &S,
				std::vector<glm::vec3> &gradS,std::vector<float> &P,
				std::vector<glm::vec3> &gradP,std::vector<glm::vec3> &gradW,
				std::vector<glm::vec3> &pie,std::vector<float> &W,float h,int sphereCount,float epsilon,
				float c,float roh0,float alpha,float jumpN,float dt,float etaMax);
#endif
		};
	} /* end namespace deformables */
} /* end namespace megamol */
#endif
