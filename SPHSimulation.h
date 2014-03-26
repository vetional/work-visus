/*
* SPHSimulation.h
*
* Copyright (C) 2011 by VISUS (Universitaet Stuttgart)
* Alle Rechte vorbehalten.
*/

#ifndef MEGAMOLCORE_SPHSimulation_H_INCLUDED
#define MEGAMOLCORE_SPHSimulation_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "CalleeSlot.h"
#include "view/AnimDataModule.h"
#include "vislib/memutils.h"
#include "vector"

//#define ONGPU 1
//#define MODELFROMFILE 1
#define DYNAMICTIME 1
#define COLLISION 1
//#define eg 1

//std includes
#include "cmath"
#include "iostream"
#include "vector"
#include <math.h> 

//using bullet physics for collision detection 
#include"btBulletCollisionCommon.h"
#include "btBulletDynamicsCommon.h"

#include "CallTriMeshData.h"


#ifdef ONGPU
#define GLM_FORCE_CUDA
#include "InitGPU.cuh"
#else
// Include GLM
#include "glm-0.9.4.6\glm\glm\glm.hpp"
#include "glm-0.9.4.6\glm/glm/gtc/matrix_transform.hpp"
#include "glm-0.9.4.6\glm/glm/gtc/quaternion.hpp"
#include "glm-0.9.4.6\glm/glm/gtx/quaternion.hpp"
#include "glm-0.9.4.6\glm/glm/gtx/euler_angles.hpp"
#include "glm-0.9.4.6\glm/glm/gtx/norm.hpp"
#include "glm-0.9.4.6/glm/glm\gtc\type_ptr.hpp"
#endif

//thrust includes
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/generate.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/version.h>

//#include "constants.h"

using namespace megamol;
using namespace megamol::core;
//using namespace constant;
namespace megamol {
	namespace deformables {


		class sideofPlane{
		public:
			int* a;
		};
		class cutPlane{
		public:
			std::vector<glm::vec3> points;
			float a;
			float b;
			float c;
			float d;
			glm::vec3 normal;
		};
		class Segment{
		public:
			glm::vec3 P0;
			glm::vec3 P1;
		};
		class treeNode{
		public:
			treeNode *leftChild;
			treeNode *rightChild;
			glm::vec4 location;
			std::vector <glm::vec4> points;
			float radius;
			int index;
		};

		/**
		* Test data source module providing generated spheres data
		*/
		class SPHSimulation : public view::AnimDataModule {
		public:

			/**
			* Answer the name of this module.
			*
			* @return The name of this module.
			*/
			static const char *ClassName(void) {
				return "SPHSimulation";
			}

			/**
			* Answer a human readable description of this module.
			*
			* @return A human readable description of this module.
			*/
			static const char *Description(void) {
				return "Provides an implemantation of Suspended particle hydrodynamics  Particle based viscoplastic fluid / solid simulation  Afonso Pavia et. al";
			}

			/**
			* Answers whether this module is available on the current system.
			*
			* @return 'true' if the module is available, 'false' otherwise.
			*/
			static bool IsAvailable(void) {
				return true;
			}

			/** Ctor. */
			SPHSimulation(void);

			/** Dtor. */
			virtual ~SPHSimulation(void);
			float *vertices;
			long time;
			std::vector<glm::vec3> triList;
			static const float globalRadius;
			std::vector<float> m;
			std::vector<float> roh;
			std::vector<float> proh;

			std::vector<float> droh;
			std::vector<glm::mat3> D;
			std::vector<float> d;
			std::vector<float>eta;


			std::vector<glm::vec3> pos;
			std::vector<glm::vec3> vel;
			std::vector<glm::mat3> gradV;
			std::vector<glm::vec3> acc;
			std::vector<glm::vec3> pacc;
			std::vector<glm::mat3> S;
			std::vector<glm::vec3> gradS;
			std::vector<float> P;
			std::vector<glm::vec3> gradP;
			std::vector<glm::vec3> gradW;
			std::vector<glm::vec3> pie;
			std::vector<float> W;
#ifdef eg
			float* h_m;
			float* h_roh;
			float* h_proh;
			float* h_droh;
			glm::mat3* h_D;
			float* h_d;
			float*h_eta;
			glm::vec3* h_pos;
			glm::vec3 * h_vel;
			glm::mat3 * h_gradV;
			glm::vec3 * h_acc;
			glm::vec3 * h_pacc;
			glm::mat3 * h_S;
			glm::vec3 * h_gradS;
			float * h_P;
			glm::vec3 * h_gradP;
			glm::vec3 * h_pie;
#endif
#ifdef ONGPU
			GPUimplementation gpu;

#endif
			bool SPHSimulation::loadOBJ(const char * path, std::vector<glm::vec3> & out_vertices);

			void SPHSimulation::initData(std::vector<float> &m,std::vector<float> &roh,std::vector<float> &proh,
				std::vector<float> &droh,std::vector<glm::mat3> &D,std::vector<float> &d,
				std::vector<float>&eta,std::vector<glm::vec3> &pos,std::vector<glm::vec3> &vel,
				std::vector<glm::mat3> &gradV,std::vector<glm::vec3> &acc,
				std::vector<glm::vec3> &pacc,std::vector<glm::mat3> &S,
				std::vector<glm::vec3> &gradS,std::vector<float> &P,
				std::vector<glm::vec3> &gradP,std::vector<glm::vec3> &gradW,
				std::vector<glm::vec3> &pie,std::vector<float> &W);

			void SPHSimulation::updategradV(std::vector<float> &m,std::vector<float> &roh,
				std::vector<glm::vec3> &vel,std::vector<glm::mat3> &gradV,
				std::vector<float> &gradW, int index,std::vector<int> &nxi);

			void SPHSimulation::calcGradW(std::vector <float>&W,
				std::vector <glm::vec3>&gradW,std::vector<glm::vec3> &pos,
				int i,std::vector<int> &nxi);

			float SPHSimulation::w(float q);

			void SPHSimulation::findNeighbours(int index, std::vector<glm::vec3> &pos,
				std::vector<int> &nxi,float h);

			void SPHSimulation::updategradV(std::vector<float> &m,std::vector<float> &roh,
				std::vector<glm::vec3> &vel,std::vector<glm::mat3> &gradV,
				std::vector<glm::vec3> &gradW, int index,std::vector<int> &nxi);

			float SPHSimulation::compute_trace(glm::mat3 D1);

			void SPHSimulation::updateD(std::vector<glm::mat3> &gradV,std::vector<glm::mat3> &D,int index);
			void SPHSimulation::updateDDE(std::vector<glm:: mat3> &D,std::vector<float> &droh,int index);//see equaqtion 10

			void SPHSimulation::compute_roh(std::vector<float>&roh,std::vector<float>&proh,std::vector<float>&m,
				std::vector<float> &W, int index,
				std::vector<int> &nxi);

			//void SPHSimulation::compute_d(std::vector <glm::mat3>D, std::vector<float> &d);

			void SPHSimulation::updateS(std::vector<float> &eta,std::vector<glm::mat3>D,
				std::vector<glm::mat3>&S,int index);

			void SPHSimulation::updateP(std::vector<float> &P,float &roh0,std::vector<float> &roh,float &c, int index);

			void SPHSimulation::updateEta(std::vector<float>&eta,std::vector<glm::mat3>&D,float n,int index);//see equaqtion 4

			void SPHSimulation::compute_gradP(std::vector<float>&m,std::vector<float>&roh,
				std::vector<float>&P,std::vector<glm::vec3>&gradW,
				std::vector<glm::vec3>&gradP,int index,std::vector<int> &nxi);

			void SPHSimulation::compute_gradS(std::vector<float>&roh,std::vector<glm::mat3>&S,
				std::vector<glm::vec3> &gradW, int index,std::vector<float> &m,
				std::vector<glm::vec3> &gradS,std::vector<int> &nxi);

			void SPHSimulation::compute_pi(std::vector<float>&m,std::vector<float>&roh,std::vector<glm::vec3>&pie,float h,
				std::vector<glm::vec3>&gradW,int index,std::vector<int> &nxi,std::vector<glm::vec3>&vel,
				std::vector<glm::vec3>&pos,float alpha,float c);

			void SPHSimulation::updateAcc(std::vector<glm::vec3>&gradP,std::vector<glm::vec3> &gradS,
				std::vector<glm::vec3>&pie,std::vector<glm::vec3> &acc,
				std::vector<glm::vec3> &pacc);//see equation 2,12

			void SPHSimulation::updateLFS(std::vector<glm::vec3> &vel,std::vector<float> &roh,std::vector<float> &proh,
				std::vector<glm::vec3> &acc,std::vector<glm::vec3> &pos,float dt,
				std::vector<glm::vec3> &pacc);//see leap frog scheme

			void SPHSimulation::correctV(std::vector<glm::vec3> &vel,std::vector<float> &m,
				std::vector<float> &W,std::vector<float> &roh,
				std::vector<glm::vec3> &pos);//see equation 11		

			/*void myCdaCallback(btBroadphasePair& collisionPair,
			btCollisionDispatcher& dispatcher, btDispatcherInfo& dispatchInfo);
			*/
			void SPHSimulation::eventualCollision(std::vector<glm::vec3> &pos,std::vector<glm::vec3> &vel);


			void SPHSimulation::updateDT(float &dt,float vmax,float etaMax,float c, float h);//see equaqtion 13

			void SPHSimulation::cudaTest();

			void SPHSimulation::computePlane(cutPlane &pl);

			void SPHSimulation::addCutPlane(std::vector<glm::vec3> &pos,std::vector<sideofPlane> &sc,
				std::vector<glm::vec3> &pt);

			bool SPHSimulation::intersect3D_RayTriangle( Segment &R, glm::vec3* T);

			bool SPHSimulation::sphereIntersectionTriangle3D(glm::vec3* T,glm::vec3 P,
				float radius);

			void SPHSimulation::broadPhase(treeNode* node, std::vector<treeNode*> &broadlist,
				std::vector<glm::vec3> triList);

			struct treeNode* SPHSimulation::computeKdTree(std::vector<glm::vec4> &points,
				treeNode *node,int depth);

			void SPHSimulation::sortPointList(std::vector <glm::vec4> &points, int axis);

			struct treeNode* SPHSimulation::NewNode();

			float SPHSimulation::boundingRadius(glm::vec4 median,std::vector <glm::vec4> &points);

			bool SPHSimulation::getTriExtentCallback(Call& caller) ;
			bool SPHSimulation::getTrigetData(Call& caller) ;
			long planeChange;
		protected:

			/**
			* Creates a frame to be used in the frame cache. This method will be
			* called from within 'initFrameCache'.
			*
			* @return The newly created frame object.
			*/
			virtual Frame* constructFrame(void) const;

			/**
			* Implementation of 'Create'.
			*
			* @return 'true' on success, 'false' otherwise.
			*/
			virtual bool create(void);

			/**
			* Loads one frame of the data set into the given 'frame' object. This
			* method may be invoked from another thread. You must take 
			* precausions in case you need synchronised access to shared 
			* ressources.
			*
			* @param frame The frame to be loaded.
			* @param idx The index of the frame to be loaded.
			*/
			virtual void loadFrame(Frame *frame, unsigned int idx);

			/**
			* Implementation of 'Release'.
			*/
			virtual void release(void);

			/**
			* Gets the data from the source.
			*
			* @param caller The calling call.
			*
			* @return 'true' on success, 'false' on failure.
			*/
			bool getDataCallback(Call& caller);

			/**
			* Gets the data from the source.
			*
			* @param caller The calling call.
			*
			* @return 'true' on success, 'false' on failure.
			*/
			bool getExtentCallback(Call& caller);

		private:

			/** Number of frames to be generated */
			static const unsigned int frameCount;

			/** Number of spheres to be generated */
			static  unsigned int sphereCount;
			// global radius of particles;

			/**
			* Class storing data of a single frame
			*/
			class Frame : public view::AnimDataModule::Frame{
			public:

				/**
				* Ctor
				*
				* @param owner The owning module
				*/
				Frame(SPHSimulation& owner) : view::AnimDataModule::Frame(owner), data(NULL) {
					// intentionally empty
				}

				/**
				* Dtor
				*/
				virtual ~Frame(void) {
					ARY_SAFE_DELETE(this->data);
				}

				/**
				* Sets the frame number
				*
				* @param n The frame number
				*/
				void SetFrameNumber(unsigned int n) {
					this->frame = n;
				}

				/** The particle data */
				float *data;

			};

			/** The slot for requesting data */
			CalleeSlot getDataSlot;

			/** The slot for requesting data */
			CalleeSlot getTriDataSlot;

			trisoup::CallTriMeshData::Mesh meshClipPlane;
		};




	} /* end namespace deformables */
} /* end namespace megamol */

#endif /* MEGAMOLCORE_SPHSimulation_H_INCLUDED */
