/*
* SPHSimulation.cpp
*
* Copyright (C) 2011 by Universitaet Stuttgart (VISUS). 
* Alle Rechte vorbehalten.
*/

#include "stdafx.h"
#define _USE_MATH_DEFINES

#include "SPHSimulation.h"
#include "moldyn/MultiParticleDataCall.h"
#include "vislib/Vector.h"
#include "vislib/ShallowVector.h"
#include "vislib/Quaternion.h"


using namespace megamol;
using namespace megamol::core;

namespace megamol
{
	namespace deformables
	{

		/*
		* SPHSimulation::frameCount
		*/
		const unsigned int SPHSimulation::frameCount = 100;
#ifdef ONGPU
		const float SPHSimulation::globalRadius=.01;
#else
		const float SPHSimulation::globalRadius=.03;
#endif

		/*
		* SPHSimulation::sphereCount
		*/


#ifdef ONGPU
		unsigned int SPHSimulation::sphereCount = 200;
#else
		unsigned int SPHSimulation::sphereCount = 150;
#endif

		/*
		* SPHSimulation::SPHSimulation
		*/

		SPHSimulation::SPHSimulation(void) : view::AnimDataModule(), getDataSlot("getData", "Gets the data from the data source"),getTriDataSlot("getTriData", "Gets the data from the data source") {
			SPHSimulation::initData(m,roh,proh,droh,D,d,eta,pos,vel,gradV,acc,pacc,S,gradS,P,gradP,gradW,pie,W);
#ifdef ONGPU
			cudaTest();
#endif			
			SPHSimulation::planeChange=0;
			this->getDataSlot.SetCallback(moldyn::MultiParticleDataCall::ClassName(), moldyn::MultiParticleDataCall::FunctionName(0), &SPHSimulation::getDataCallback);
			this->getDataSlot.SetCallback(moldyn::MultiParticleDataCall::ClassName(), moldyn::MultiParticleDataCall::FunctionName(1), &SPHSimulation::getExtentCallback);
			this->MakeSlotAvailable(&this->getDataSlot);
			this->getTriDataSlot.SetCallback(trisoup::CallTriMeshData::ClassName(), trisoup::CallTriMeshData::FunctionName(0), &SPHSimulation::getTrigetData);
			this->getTriDataSlot.SetCallback(trisoup::CallTriMeshData::ClassName(), trisoup::CallTriMeshData::FunctionName(1), &SPHSimulation::getTriExtentCallback);
			this->MakeSlotAvailable(&this->getTriDataSlot);


			float *vertices = new float[18];

			vertices[0] = 0.0f; vertices[1] = 0.0f; vertices[2] = 0.0f;
			vertices[3] = 0.0f; vertices[4] = 2.50f; vertices[5] = 0.0f;
			vertices[6] = 0.0f; vertices[7] = 0.0f; vertices[8] = 02.50f;
			vertices[9] = 0.0f; vertices[10] = 0.0f; vertices[11] = 0.0f;
			vertices[12] = 0.0f; vertices[13] = -1.50f; vertices[14] = 0.0f;
			vertices[15] = 0.0f; vertices[16] = 0.0f; vertices[17] = -02.50f;

			for (int i = 0; i < 18; i+=3)
			{
				triList.push_back(glm::vec3 (vertices[i],vertices[i+1],vertices[i+2])); 
			}

			meshClipPlane.SetVertexData<float*, float*, float*, float*>(6, vertices, NULL, NULL, NULL, true);

		}


		/*
		* SPHSimulation::~SPHSimulation
		*/
		SPHSimulation::~SPHSimulation(void) {
			this->Release();
		}


		/*
		* SPHSimulation::constructFrame
		*/
		view::AnimDataModule::Frame* SPHSimulation::constructFrame(void) const {
			return new Frame(const_cast<SPHSimulation&>(*this));
		}


		/*
		* BezierDataSource::create
		*/
		bool SPHSimulation::create(void) {
			view::AnimDataModule::setFrameCount(SPHSimulation::frameCount);
			view::AnimDataModule::initFrameCache(SPHSimulation::frameCount);
			return true;
		}


		/*
		* SPHSimulation::loadFrame
		*/
		void SPHSimulation::loadFrame(view::AnimDataModule::Frame *frame, unsigned int idx) {
			Frame *frm = dynamic_cast<Frame *>(frame);
			if (frm == NULL) return;
			frm->SetFrameNumber(idx);
			frm->data = new float[7 * SPHSimulation::sphereCount];
			for (unsigned int i = 0; i < SPHSimulation::sphereCount; i++) {
				vislib::math::ShallowVector<float, 3> pos(&frm->data[i * 7]);
				::srand(i); // stablize values for particles
				float &r = frm->data[i * 7 + 3];
				float &cr = frm->data[i * 7 + 4];
				float &cg = frm->data[i * 7 + 5];
				float &cb = frm->data[i * 7 + 6];
				vislib::math::Vector<float, 3> X(static_cast<float>((::rand() % 2) * 2 - 1), 0.0f, 0.0f);
				vislib::math::Vector<float, 3> Y(0.0f, static_cast<float>((::rand() % 2) * 2 - 1), 0.0f);
				vislib::math::Vector<float, 3> Z(
					static_cast<float>(1000 - ::rand() % 2001) * 0.001f,
					static_cast<float>(1000 - ::rand() % 2001) * 0.001f,
					static_cast<float>(1000 - ::rand() % 2001) * 0.001f);
				switch (::rand() % 6) {
				case 0: Z.SetX(1.0f); break;
				case 1: Z.SetX(-1.0f); break;
				case 2: Z.SetY(1.0f); break;
				case 3: Z.SetY(-1.0f); break;
				case 4: Z.SetZ(1.0f); break;
				case 5: Z.SetZ(-1.0f); break;
				}
				Z.Normalise();
				vislib::math::Quaternion<float> rot(static_cast<float>((::rand() % 2) * 2 - 1) * static_cast<float>(M_PI) * static_cast<float>(::rand() % 2000) * 0.001f, Z);
				float dist = static_cast<float>(::rand() % 1001) * 0.001f;
				dist = ::pow(dist, 0.333f) * 0.9f;
				float a = (static_cast<float>(2 * idx) / static_cast<float>(SPHSimulation::frameCount)) * static_cast<float>(M_PI);
				X = rot * X;
				Y = rot * Y;

				X *= sin(a) * dist;
				Y *= cos(a) * dist;
				pos = X;
				pos += Y;

				r = 0.05f + static_cast<float>(::rand() % 501) * 0.0001f;

				Z.Set(
					static_cast<float>(1000 - ::rand() % 2001) * 0.001f,
					static_cast<float>(1000 - ::rand() % 2001) * 0.001f,
					static_cast<float>(1000 - ::rand() % 2001) * 0.001f);
				switch (::rand() % 6) {
				case 0: Z.SetX(1.0f); break;
				case 1: Z.SetX(-1.0f); break;
				case 2: Z.SetY(1.0f); break;
				case 3: Z.SetY(-1.0f); break;
				case 4: Z.SetZ(1.0f); break;
				case 5: Z.SetZ(-1.0f); break;
				}
				Z.Normalise();

				cr = vislib::math::Abs(Z.X());
				cg = vislib::math::Abs(Z.Y());
				cb = vislib::math::Abs(Z.Z());

			}

		}


		/*
		* SPHSimulation::release
		*/
		void SPHSimulation::release(void) {
			// intentionally empty
		}


		/*
		* SPHSimulation::getDataCallback
		*/
		//--------------------------------------------------------------------------------------------------------------------------
		// Global DATA Store for Evaluation Function
		//--------------------------------------------------------------------------------------------------------------------------
		long time=0;
		float h=20*SPHSimulation::globalRadius;// kernel radius
		float dt=0.0000001f;//time step
		float vmax=150;//for dt update
		float etaMax=1600;//for dt update viscosity of the materail is represented by eta.
		float alpha=1.5;//bulk viscosity
		float epsilon=0.30f;//for velocity correction due to xsph value=[0,1]
		float roh0=8.0;//reference density -- the system will calculate pressure with this reference and if not set correctly it might greatly deform the system.
		float c=1030;// speed of sound in the medium in m/sec
		float n=0.5f;//refer equations in paper for eta
		float jumpN=12.5f;//jump number
		glm::vec3 velocity(-150.0,-10,0.);//initial velocity of each particle
		float commonMass=1;// common mass of each point in the material
		float densityInit=100;//initial density vector initializied to this value.


		int numCut=0;
		std::vector<glm::vec3> triList;
		std::vector<sideofPlane> sc;
		std::vector<float> m;
		std::vector<float> roh;//stores density values for each particle
		std::vector<float> proh;//stores density values for each particle in last timestep
		std::vector<float> droh;//stores derrivative of density values for each particle
		std::vector<glm::mat3> D;//velocity gradient vector
		std::vector<float> d;//needed to calculate eta
		std::vector<float>eta;//viscosity of each particle
		std::vector<glm::vec3> pos;// position of each particle
		std::vector<glm::vec3> vel;//velocity of each particle
		std::vector<glm::mat3> gradV;//gradient of velocity for each particle
		std::vector<glm::vec3> acc;// accleration of each particle
		std::vector<glm::vec3> pacc;//accleration of each particle in last timestep
		std::vector<glm::mat3> S;//
		std::vector<glm::vec3> gradS;//
		std::vector<float> P;//pressure at each particle
		std::vector<glm::vec3> gradP;//gradient of pressure at each particle
		std::vector<glm::vec3> gradW;//gradient of kernel for each of the neighbour of the particle
		std::vector<glm::vec3> pie;
		std::vector<float> W;//kernel value for each neighbour of the particle
		glm::vec3 z(0.0f,0.0f,0.0f);
		glm::mat3 zero(z,z,z);
		float ze=0.0f;
		btVector3 fallInertia(0,0,0);
#ifdef ONGPU
		GPUimplementation gpu;

#endif
		//--------------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------------

		void SPHSimulation::initData(std::vector<float> &m,std::vector<float> &roh,std::vector<float> &proh,
			std::vector<float> &droh,std::vector<glm::mat3> &D,std::vector<float> &d,
			std::vector<float>&eta,std::vector<glm::vec3> &pos,
			std::vector<glm::vec3> &vel,
			std::vector<glm::mat3> &gradV,std::vector<glm::vec3> &acc,
			std::vector<glm::vec3> &pacc,std::vector<glm::mat3> &S,
			std::vector<glm::vec3> &gradS,std::vector<float> &P,
			std::vector<glm::vec3> &gradP,std::vector<glm::vec3> &gradW,
			std::vector<glm::vec3> &pie,std::vector<float> &W){



				for(int i = 0;i<SPHSimulation::sphereCount;i++){
					m.push_back(commonMass);
					roh.push_back(densityInit);
					proh.push_back(densityInit);
					droh.push_back(densityInit);
					D.push_back(zero);
					d.push_back(ze);

					vel.push_back(velocity);
					gradV.push_back(zero);
					acc.push_back(z);
					pacc.push_back(z);
					S.push_back(zero);
					gradS.push_back(z);
					P.push_back(ze);
					gradP.push_back(z);
					pie.push_back(z);

					eta.push_back(ze);
#ifndef MODELFROMFILE
					pos.push_back(z);
#endif
				}
				float a1=0.3;float b1=.15;float c1=0.;
				int dc=0;
#ifndef ONGPU
				int side=5;
#else
				int side=5;
#endif
				while(dc<SPHSimulation::sphereCount){



					for (int ja = 0; ja < side; ja++)
					{
						for (int ka = 0; ka < side; ka++)
						{
							if(dc<SPHSimulation::sphereCount){
								pos[dc].x=a1;
								pos[dc].y=b1;
								pos[dc].z=c1;
								dc++;
							}
							a1+=globalRadius*8.5f;
						}
						a1=0.3;
						b1+=globalRadius*6.5f;

					}c1+=globalRadius*8.5f;b1=.15;

				}
#ifdef MODELFROMFILE
				//	bool res_obj = loadOBJ("H:\\MEGAmol\\deformables\\trunk\\models\\bunny01.obj", pos);
				bool res_obj = loadOBJ("H:\\MEGAmol\\deformables\\trunk\\models\\sphere_lowres.obj", pos);
#endif
				btCollisionShape* fallShape = new btSphereShape(globalRadius);
				fallShape->calculateLocalInertia(commonMass,fallInertia);
		}
		bool SPHSimulation::getDataCallback(Call& caller) {
			moldyn::MultiParticleDataCall *mpdc = dynamic_cast<moldyn::MultiParticleDataCall *>(&caller);
			if (mpdc == NULL) return false;

			view::AnimDataModule::Frame *f = this->requestLockedFrame(mpdc->FrameID());
			if (f == NULL) return false;
			f->Unlock(); // because I know that this data source is simple enough that no locking is required
			Frame *frm = dynamic_cast<Frame*>(f);
			if (frm == NULL) return false;



			//--------------------------------------------------------------------------------------------------------------------------
			// update-data using sph formulation.
			//--------------------------------------------------------------------------------------------------------------------------

#ifndef ONGPU
			//See equations from Paper : Particle based viscoplastic fluid / solid simulation  Afonso Pavia et. al 
			for(int i = 0;i<m.size();i++){

				std::vector<int> nxi;
				SPHSimulation::findNeighbours(i,pos,nxi,h);
				//nxi containts the list of neighbours of i

				SPHSimulation::calcGradW(W,gradW,pos,i,nxi);// calculate w and gradW as they are required for each of next fn calls
				// Turns out we don't need to calcultate w for gradW
				//but we need w for roh
				SPHSimulation::compute_roh(roh,proh,m,W,i,nxi);
				SPHSimulation::updategradV(m,roh,vel,gradV,gradW,i,nxi);//see equqtion 7

				SPHSimulation::updateD(gradV,D,i); //see equaqtion 8

				SPHSimulation::updateDDE(D,droh,i);//see equaqtion 10

				//compute roh 
				SPHSimulation::compute_roh(roh,proh,m,W,i,nxi);
				SPHSimulation::updateP(P,roh0,roh,c,i);//see equaqtion 6

				//SPHSimulation::compute_d(D,d); // extra function to compute trace of D
				SPHSimulation::updateEta(eta,D,n,i);//see equaqtion 4

				SPHSimulation::updateS(eta,D,S,i);

				//for artificial viscosity equation 12
				SPHSimulation::compute_pi(m,roh,pie,h,gradW,i,nxi,vel,pos,alpha,c);

				// to update acceleration we need to calculate grad p and grad s 
				SPHSimulation::compute_gradS(roh,S,gradW,i,m,gradS,nxi);
				SPHSimulation::compute_gradP(m,roh,P,gradW,gradP,i,nxi);
			}

			SPHSimulation::updateAcc(gradP,gradS,pie,acc,pacc);//see equaqtion 2,12

			SPHSimulation::updateLFS(vel,roh,proh,acc,pos,dt,pacc);//see leap frog scheme
			SPHSimulation::correctV(vel,m,W,roh,pos);//see equaqtion 11
#endif
#ifdef ONGPU

			gpu.theBigLoop(m,roh,proh,droh,D,d,eta,pos,vel,gradV,acc,pacc,S,gradS,P,
				gradP,gradW,pie,W, h,sphereCount,epsilon,c,roh0,alpha,jumpN,dt,etaMax);
#endif
#ifdef COLLISION
			SPHSimulation::eventualCollision(pos,vel);
#endif
#ifdef DYNAMICTIME
			SPHSimulation::updateDT(dt,vmax,etaMax,c,h);//see equaqtion 13
#endif
			time+=dt;



			//--------------------------------------------------------------------------------------------------------------------------

			//--------------------------------------------------------------------------------------------------------------------------

			if(pos.size() > 0)
			{
				mpdc->SetFrameID(f->FrameNumber());
				mpdc->SetDataHash(1);
				mpdc->SetExtent(SPHSimulation::frameCount,
					-1.0f, -1.0f, -1.0f,
					1.0f, 1.0f, 1.0f);
				mpdc->SetParticleListCount(1);
				mpdc->AccessParticles(0).SetGlobalRadius(SPHSimulation::globalRadius);
				mpdc->AccessParticles(0).SetCount(SPHSimulation::sphereCount);
				mpdc->AccessParticles(0).SetVertexData(moldyn::SimpleSphericalParticles::VERTDATA_FLOAT_XYZ, pos.data(), sizeof(float) * 3);
				mpdc->AccessParticles(0).SetColourData(moldyn::SimpleSphericalParticles::COLDATA_FLOAT_RGB, frm->data + 4, sizeof(float) * 7);
				mpdc->SetUnlocker(NULL);
			}
			return true;
		}


		bool SPHSimulation::getTriExtentCallback(Call& caller) {

			trisoup::CallTriMeshData *ctmd = dynamic_cast<trisoup::CallTriMeshData *>(&caller);
			if (ctmd == NULL) return false;

			ctmd->SetDataHash(1);
			ctmd->SetExtent(SPHSimulation::frameCount,
				-1.0f, -1.0f, -1.0f,
				1.0f, 1.0f, 1.0f);
			return true;
		}
		bool SPHSimulation::getTrigetData(Call& caller) {

			trisoup::CallTriMeshData *ctmd = dynamic_cast<trisoup::CallTriMeshData *>(&caller);
			if (ctmd == NULL) return false;


			ctmd->SetDataHash(SPHSimulation::planeChange++);
			ctmd->SetObjects(1, &meshClipPlane);
			ctmd->SetUnlocker(NULL);

			return true;
		}
		/*
		* SPHSimulation::getExtentCallback
		*/
		bool SPHSimulation::getExtentCallback(Call& caller) {
			moldyn::MultiParticleDataCall *mpdc = dynamic_cast<moldyn::MultiParticleDataCall *>(&caller);
			if (mpdc == NULL) return false;

			mpdc->SetDataHash(1);
			mpdc->SetExtent(SPHSimulation::frameCount,
				-1.0f, -1.0f, -1.0f,
				1.0f, 1.0f, 1.0f);

			return true;
		}

		float SPHSimulation::compute_trace(glm::mat3 D1){
			return D1[0][0]+D1[1][1]+D1[2][2];
		}

		//calculates w but is not needed for calculating gradW 
		//update need w for roh calculation.
		float SPHSimulation::w(float q)
		{
			float res=0;
			if(q<=2||q>=0)
			{
				res=0.66+((9/8)*q*q)+((19/24)*q*q*q)-((5/32)*q*q*q*q);
			}
			return res;
		}

		void  SPHSimulation::computePlane(cutPlane &pl)
		{


			float vx = (pl.points[1].x - pl.points[2].x);
			float vy = (pl.points[1].y - pl.points[2].y);
			float vz = (pl.points[1].z - pl.points[2].z);

			float wx = (pl.points[0].x - pl.points[1].x);
			float wy = (pl.points[0].y - pl.points[1].y);
			float wz = (pl.points[0].z - pl.points[1].z);

			float vw_x = vy * wz - vz * wy;
			float vw_y = vz * wx - vx * wz;
			float vw_z = vx * wy - vy * wx;

			float mag = sqrtf((vw_x * vw_x) + (vw_y * vw_y) + (vw_z * vw_z));

			if ( mag > 0.000001f )
			{

				mag = 1.0f/mag; // compute the recipricol distance

				pl.a = vw_x * mag;
				pl.b = vw_y * mag;
				pl.c = vw_z * mag;
				pl.d = 0.0f - ((pl.a*pl.points[0].x)+(pl.b*pl.points[0].y)+(pl.c*pl.points[0].z));

			}

		}
		void SPHSimulation::addCutPlane(std::vector<glm::vec3> &pos,std::vector<sideofPlane> &sc,std::vector<glm::vec3> &pt){

			cutPlane plane;

			for (int i = 0; i < 4; i++)
			{
				plane.points.push_back(pt[i]);
			}

			computePlane(plane);

			for (int i = 0; i < SPHSimulation::sphereCount; i++)
			{
				//a*x+b*y+c*z+d = 0;
				sideofPlane pl;
				sc.push_back(pl);

				if(glm::dot( glm::vec4(plane.a,plane.b,plane.c,plane.d), glm::vec4(pos[i].x,pos[i].y,pos[i].z,1) )>0 ){
					sc[i].a[numCut]=numCut;
				}
				else{
					sc[i].a[numCut]=-numCut;
				}
			}

		}
		void SPHSimulation::findNeighbours(int index, std::vector<glm::vec3> &pos, 
			std::vector<int> &nxi,float h){

				for(int i =0; i <pos.size();i++){

					if(glm::distance(pos[index],pos[i])<h){

						if(index!=i){

							int flag=0;

							for (int j = 0; j < numCut; j++)
							{
								if(sc[i].a[j]!=sc[index].a[j]){
									flag=1;
									break;
								}
							}
							if(flag==0)
								nxi.push_back(i);

						}
					}

				}

		}
		// the pre requisite to call this function is the nearest neighbours should 
		// have been already computed.

		void SPHSimulation::calcGradW(std::vector <float>&W,
			std::vector <glm::vec3>&gradW,std::vector<glm::vec3> &pos,
			int i,std::vector<int> &nxi){


				//using w(q) .66+9/8*q^2+19/24*q^3-5/32*q^4
				// W(x-xj,h)= 315/208*3.14*h^3 w(mod(x-xj)/h)

				W.clear();
				gradW.clear();
				for(int j=0;j<nxi.size();j++)
				{
					W.push_back(ze);
					gradW.push_back(z);
					W[j]=((315/208)*3.14*h*h*h)* w(glm::distance(pos[i],pos[nxi[j]])/h);
					float dXiXj=glm::distance(pos[i],pos[nxi[j]]);
					float multiplier=-(9/(4*h))+(19/(8*h*h))*dXiXj+(5/(8*h*h*h))*dXiXj;

					gradW[j].x= multiplier*((float)pos[i].x-(float)pos[nxi[j]].x);
					gradW[j].y= multiplier*(pos[i].y-pos[nxi[j]].y);
					gradW[j].z= multiplier*(pos[i].z-pos[nxi[j]].z);
				}
		}


		void SPHSimulation::updategradV(std::vector<float> &m,std::vector<float> &roh,
			std::vector<glm::vec3> &vel,std::vector<glm::mat3> &gradV,
			std::vector<glm::vec3> &gradW, int index,std::vector<int> &nxi)
		{

			for(int j=0;j<nxi.size();j++)
			{

				glm::vec3 v=(vel[nxi[j]]-vel[index]);
				gradV[index]+=((m[nxi[j]]/roh[nxi[j]])*glm::outerProduct(v,gradW[j]));

			}

		}
		void SPHSimulation::updateD(std::vector<glm::mat3> &gradV,std::vector<glm::mat3> &D,int index){

			D[index]=gradV[index]+glm::transpose(gradV[index]);
		}
		void SPHSimulation::updateDDE(std::vector<glm:: mat3> &D,std::vector<float> &droh,int index){


			glm::mat3 D1= 0.5*D[index];
			droh[index] = D1[0][0]+D1[1][1]+D1[2][2];

		}//see equaqtion 10
		void SPHSimulation::compute_roh(std::vector<float>&roh,std::vector<float>&proh,std::vector<float>&m,
			std::vector<float> &W, int index,std::vector<int> &nxi){

				float sum=0;
				for(int j=0;j<nxi.size();j++)
				{
					sum+=m[nxi[j]]*W[j];
				}
				proh[index]=roh[index];
				roh[index]=sum;
		}

		//|Updateds the pressure
		//Todo; : check the memeory allcation of the pressure term before the function call.
		void SPHSimulation::updateP(std::vector<float> &P,float &roh0,std::vector<float> &roh,
			float &c, int index){

				P[index]=c*c*(roh[index]-roh0);

		}//see equaqtion 6

		void SPHSimulation::updateEta(std::vector<float>&eta,std::vector<glm::mat3>&D,float n,int index) //see equaqtion 4
		{

			glm::mat3 D1= D[index];

			//compute the d from the current rate of deformation tensor 
			float d11=std::sqrt(0.5*(D1[0][0]+D1[1][1]+D1[2][2])*(D1[0][0]+D1[1][1]+D1[2][2]));
			//the whole system crashes as d11 is zero initially 
			//due to same velocity of each particle at startup.
			//therefore counter measure taken.(has to be verified though)
			if(d11==0)
				d11=1;

			float exp=std::exp(-(d11*(jumpN+1)));
			// since we have fixed n to be .5 other wise use the commented version.
			float temp=(1-exp)*((1/d11)+(1/d11));
			if(temp<etaMax)
				eta[index]=temp;
			else
				eta[index]=etaMax;
			//eta[index]=(1-exp)*(std::pow(d11,n-1)+(1/d11));

		}

		void SPHSimulation::updateS(std::vector<float> &eta,std::vector<glm::mat3>D,
			std::vector<glm::mat3>&S,int index){


				S[index]=eta[index]*D[index];

		}
		void SPHSimulation::compute_pi(std::vector<float>&m,std::vector<float>&roh,std::vector<glm::vec3>&pie,float h,
			std::vector<glm::vec3>&gradW,int index,std::vector<int> &nxi,std::vector<glm::vec3>&vel,
			std::vector<glm::vec3>&pos,float alpha,float c){

				glm::vec3 sum(0.0f,0.0f,0.0f);
				for(int j=0;j<nxi.size();j++)
				{

					if(glm::dot((vel[index]-vel[nxi[j]]),(pos[index]-pos[nxi[j]]))<0)
					{
						float dist=glm::distance(pos[index],pos[nxi[j]]);
						float mu=h*glm::dot((vel[index]-vel[nxi[j]]),(pos[index]-pos[nxi[j]]))/(dist*dist+0.01*h*h);
						sum+=m[nxi[j]]*((2*alpha*c*(h*mu))/(roh[index]+roh[j]))*gradW[j];
					}
					else
						sum+=0;
				}
				pie[index]=sum;

		}
		// to update acceleration we need to calculate grad p and grad s 

		void SPHSimulation::compute_gradS(std::vector<float>&roh,std::vector<glm::mat3>&S,
			std::vector<glm::vec3> &gradW, int index,std::vector<float> &m,
			std::vector<glm::vec3> &gradS,std::vector<int> &nxi){

				glm::vec3 sum(0.0f,0.0f,0.0f);
				for(int j=0;j<nxi.size();j++)
				{
					//using the formulation in paper.
					glm::mat3 sv=(m[nxi[j]]/(roh[nxi[j]]*roh[index]))*(S[nxi[j]]+S[index]);
					sum+=sv*gradW[j];

				}
				gradS[index]=sum;
		}
		void SPHSimulation::compute_gradP(std::vector<float>&m,std::vector<float>&roh,
			std::vector<float>&P,std::vector<glm::vec3>&gradW,
			std::vector<glm::vec3>&gradP,int index,std::vector<int> &nxi){

				glm::vec3 sum(0.0f,0.0f,0.0f);
				float ptemp;
				for(int j=0;j<nxi.size();j++)
				{

					ptemp=(m[nxi[j]])*((P[nxi[j]]/std::pow(roh[nxi[j]],2))+(P[index]/std::pow(roh[index],2)));
					sum+=ptemp*gradW[j];

				}
				gradP[index]=sum;
		}

		void SPHSimulation::updateAcc(std::vector<glm::vec3>&gradP,std::vector<glm::vec3> &gradS,
			std::vector<glm::vec3>&pie,std::vector<glm::vec3> &acc,std::vector<glm::vec3> &pacc){

				glm::vec3 gr(0.0f,-10.0f,0.0f);
				for(int i=0;i<gradP.size();i++){

					pacc[i]=acc[i];
					acc[i]=gradS[i]-gradP[i]+gr-pie[i];

				}

		}//see equaqtion 2,12


		//see leap frog scheme
		void SPHSimulation::updateLFS(std::vector<glm::vec3> &vel,std::vector<float> &roh,
			std::vector<float> &proh,std::vector<glm::vec3> &acc,std::vector<glm::vec3> &pos,float dt,
			std::vector<glm::vec3> &pacc){

				for(int i =0; i <pos.size();i++){

					vel[i]=vel[i]+0.5f*(pacc[i]+acc[i])*dt;
					// roh has been updated in the paper using LFS which is not reasonable to me 
					// as the density should be dependent only on the current configuration but it seems 
					// it needs some sort of help from previous stages.

					roh[i]=roh[i]+0.5f*(proh[i]+roh[i])*dt;
				}

		}
		void SPHSimulation::correctV(std::vector<glm::vec3> &vel,std::vector<float> &m,
			std::vector<float> &W,std::vector<float> &roh,
			std::vector<glm::vec3> &pos){//see equaqtion 11	

				for(int i = 0;i<m.size();i++){

					std::vector<int> nxi;
					SPHSimulation::findNeighbours(i,pos,nxi,h);
					glm::vec3 sum(0.0f,0.0f,0.0f);
					SPHSimulation::calcGradW(W,gradW,pos,i,nxi);
					for(int j=0;j<nxi.size();j++)
					{
						sum+=((2*m[nxi[j]])/(roh[i]+roh[nxi[j]]))*(vel[nxi[j]]-vel[i])*W[j];
					}		



					vel[i]=vel[i]+epsilon*sum;
					pos[i]=pos[i]+vel[i]*dt;

				}
		}


		// Assume that classes are already given for the objects:
		//    Point and Vector with
		//        coordinates {float x, y, z;}
		//        operators for:
		//            == to test  equality
		//            != to test  inequality
		//            Point   = Point ± Vector
		//            Vector =  Point - Point
		//            Vector =  Scalar * Vector    (scalar product)
		//            Vector =  Vector * Vector    (3D cross product)
		//    Line and Ray and Segment with defining  points {Point P0, P1;}
		//        (a Line is infinite, Rays and  Segments start at P0)
		//        (a Ray extends beyond P1, but a  Segment ends at P1)
		//    Plane with a point and a normal {Point V0; Vector  n;}
		//   #define SMALL_NUM   0.00000001 // anything that avoids division overflow
		////  dot product (3D) which allows vector operations in arguments
		//  #define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)
		//  #define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)  // perp product  (2D)

		//===================================================================




		// intersect3D_SegmentPlane(): find the 3D intersection of a segment and a plane
		//    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
		//    Output: *I0 = the intersect point (when it exists)
		//    Return: 0 = disjoint (no intersection)
		//            1 =  intersection in the unique point *I0
		//            2 = the  segment lies in the plane


		//int intersect3D_SegmentPlane( Segment S, cutPlane Pn, glm::vec3* I )
		//{
		//    Vector    u = S.P1 - S.P0;
		//    Vector    w = S.P0 - Pn.V0;
		//
		//    float     D = dot(Pn.n, u);
		//    float     N = -dot(Pn.n, w);
		//
		//    if (fabs(D) < SMALL_NUM) {           // segment is parallel to plane
		//        if (N == 0)                      // segment lies in plane
		//            return 2;
		//        else
		//            return 0;                    // no intersection
		//    }
		//    // they are not parallel
		//    // compute intersect param
		//    float sI = N / D;
		//    if (sI < 0 || sI > 1)
		//        return 0;                        // no intersection
		//
		//    *I = S.P0 + sI * u;                  // compute segment intersect point
		//    return 1;
		//}

		// intersect3D_RayTriangle(): find the 3D intersection of a ray with a triangle
		//    Input:  a ray R, and a triangle T
		//    Output: *I = intersection point (when it exists)
		//    Return: -1 = triangle is degenerate (a segment or point)
		//             0 =  disjoint (no intersect)
		//             1 =  intersect in unique point I1
		//             2 =  are in the same plane
		bool SPHSimulation::intersect3D_RayTriangle( Segment &R, glm::vec3* T)
		{
			glm::vec3 I ;

			glm::vec3    u, v, n;              // triangle vectors
			glm::vec3   dir, w0, w;           // ray vectors
			float     r, a, b;              // params to calc ray-plane intersect

			// get triangle edge vectors and plane normal
			u = T[1] - T[0];
			v = T[2] - T[0];
			n = u * v;              // cross product
			if (n == (glm::vec3)(0,0,0))             // triangle is degenerate
				return FALSE;                  // do not deal with this case

			dir = R.P1 - R.P0;              // ray direction vector
			w0 = R.P0 - T[0];
			a = -glm::dot(n,w0);
			b = glm::dot(n,dir);
			if (fabs(b) < .0000001) {     // ray is  parallel to triangle plane
				if (a == 0)                 // ray lies in triangle plane
					return TRUE;
				else return FALSE;              // ray disjoint from plane
			}

			// get intersect point of ray with triangle plane
			r = a / b;
			if (r < 0.0)                    // ray goes away from triangle
				return FALSE;                   // => no intersect
			// for a segment, also test if (r > 1.0) => no intersect

			I = R.P0 + r * dir;            // intersect point of ray and plane

			// is I inside T?
			float    uu, uv, vv, wu, wv, D;
			uu = glm::dot(u,u);
			uv = glm::dot(u,v);
			vv = glm::dot(v,v);
			w = I - T[0];
			wu = glm::dot(w,u);
			wv = glm::dot(w,v);
			D = uv * uv - uu * vv;

			// get and test parametric coords
			float s, t;
			s = (uv * wv - vv * wu) / D;
			if (s < 0.0 || s > 1.0)         // I is outside T
				return FALSE;
			t = (uv * wu - uu * wv) / D;
			if (t < 0.0 || (s + t) > 1.0)  // I is outside T
				return FALSE;

			return TRUE;                       // I is in T
		}
		//
		//Separating axis testing for a sphere
		//
		//From the above we conclude that a triangle has three kinds of features that are involved in testing: vertices (3), edges (3), and a face (1). But what about a sphere? Well, if you followed the justification of the test, it should be clear a sphere has only one feature, but it varies from test to test: the point on the sphere surface closest to each given triangle feature!
		//
		//If we define the problem as testing the intersection between a triangle T, defined by the vertices A, B, and C, and a sphere S, given by its center point P and a radius r, the axes we need to test are therefore the following:
		//
		//    normal to the triangle plane
		//    perpendicular to AB through P
		//    perpendicular to BC through P
		//    perpendicular to CA through P
		//    through A and P
		//    through B and P
		//    through C and P 

		bool SPHSimulation::sphereIntersectionTriangle3D(glm::vec3* T,glm::vec3 P,float radius){
			//T is the triangle and P is the center of the sphere
			glm::vec3 A = T[0] - P;
			glm::vec3 B = T[1] - P;
			glm::vec3 C = T[2] - P;
			float rr = radius*radius;
			glm::vec3 V = glm::cross(B - A, C - A);
			float d = glm::dot(A, V);
			float e = glm::dot(V, V);
			bool sep1=FALSE;
			if( d * d > rr * e)
				sep1 =TRUE;

			float aa = glm::dot(A, A);
			float ab = glm::dot(A, B);
			float ac = glm::dot(A, C);
			float bb = glm::dot(B, B);
			float bc = glm::dot(B, C);
			float cc = glm::dot(C, C);
			bool sep2 = (aa > rr) && (ab > aa) && (ac > aa);
			bool sep3 = (bb > rr) && (ab > bb) && (bc > bb);
			bool sep4 = (cc > rr) && (ac > cc) && (bc > cc);

			glm::vec3 AB = B - A;
			glm::vec3 BC = C - B;
			glm::vec3 CA = A - C;
			float d1 = ab - aa;
			float d2 = bc - bb;
			float d3 = ac - cc;
			float e1 = glm::dot(AB, AB);
			float e2 = glm::dot(BC, BC);
			float e3 = glm::dot(CA, CA);
			glm::vec3 Q1 = A * e1 - d1 * AB;
			glm::vec3 Q2 = B * e2 - d2 * BC;
			glm::vec3 Q3 = C * e3 - d3 * CA;
			glm::vec3 QC = C * e1 - Q1;
			glm::vec3 QA = A * e2 - Q2;
			glm::vec3 QB = B * e3 - Q3;
			bool sep5 = (glm::dot(Q1, Q1) > rr * e1 * e1) && (glm::dot(Q1, QC) > 0);
			bool sep6 = (glm::dot(Q2, Q2) > rr * e2 * e2) &&(glm::dot(Q2, QA) > 0);
			bool sep7 = (glm::dot(Q3, Q3) > rr * e3 * e3) && (glm::dot(Q3, QB) > 0);

			bool separated = sep1 || sep2 ||sep3 || sep4 || sep5 || sep6 || sep7;

			return !separated;
		}


		void SPHSimulation::broadPhase(treeNode* node, std::vector<treeNode*> &broadlist,
			std::vector<glm::vec3> triList)
		{
			//here the


			if(node->leftChild==NULL && node->rightChild==NULL)
			{
				bool dTN;
				for (int i = 0; i < triList.size(); i+=3)
				{

					glm::vec3 *T;
					T=(glm::vec3*) malloc (3);
					T[0]=triList[i];
					T[1]=triList[i+1];
					T[2]=triList[i+2];
					dTN=SPHSimulation::sphereIntersectionTriangle3D(T,glm::vec3(node->location.x,node->location.y,node->location.z),node->radius);

					if(dTN==TRUE)
					{

						break;
					}

				}
				if(dTN==TRUE)
				{
					broadlist.push_back(node);

				}
				std::cout<<"broadlist point inserted";
			}
			else
			{
				SPHSimulation::broadPhase(node->leftChild,broadlist,triList);
				SPHSimulation::broadPhase(node->rightChild,broadlist,triList);
			}




		}





		float SPHSimulation::boundingRadius(glm::vec4 median,std::vector <glm::vec4> &points){

			float dist=0;
			float temp;
			for(int kk=0;kk<points.size();kk++){

				temp=glm::distance(glm::vec3(median.x,median.y,median.z),glm::vec3(points[kk].x,points[kk].y,points[kk].z));

				if(dist<temp)
					dist=temp;
			}

			return dist;

		}

		struct treeNode* SPHSimulation::NewNode() {
			struct treeNode* node = new (struct treeNode);    // "new" is like "malloc"

			node->location=glm::vec4(0,0,0,0);
			node->leftChild = NULL;
			node->rightChild = NULL;
			node->radius=0;
			node->points.push_back(glm::vec4(0,0,0,0));

			return(node);
		} 

		void SPHSimulation::sortPointList(std::vector <glm::vec4> &points, int axis){



			for(int i=points.size()-1;i>0;i--){

				for(int j=0;j<i;j++){

					if(points[i][axis]<points[j][axis])
					{
						glm::vec4 temp=points[i];
						points[i]=points[j];
						points[j]=temp;
					}

				}
			}
		}

		struct treeNode* SPHSimulation::computeKdTree(std::vector<glm::vec4> &points,
			treeNode *node,int depth){



				int axis=depth%3;

				sortPointList(points,axis);
				int in=floor((float)points.size()/2);
				glm::vec4 median (points[in]);	

				if(depth==0)
				{

					node->location=median;
				}

				if(points.size()<20)
				{

					node->location=median;

					for(int fill=0;fill<points.size();fill++)
					{
						node->points.push_back(points[fill]);
					}
					node->radius=boundingRadius(median,points);
					node->leftChild=NULL;
					node->rightChild=NULL;
				}

				else
				{

					std::vector<glm::vec4> lcpoints;
					std::vector<glm::vec4> rcpoints;

					for(int i =0;i<points.size();i++)
					{

						glm::vec4 point=points[i];

						if(i<= floor((float)points.size()/2))
							lcpoints.push_back(point);
						else 
							rcpoints.push_back(point);

					}

					treeNode* emptylc=NewNode();
					treeNode* emptyrc=NewNode();

					if(depth!=0)
						node->location=median;

					node->radius=boundingRadius(median,points);

					node->leftChild=computeKdTree(lcpoints,emptylc,depth+1);
					node->rightChild=computeKdTree(rcpoints,emptyrc,depth+1);

				}
				return node;
		}

		void SPHSimulation::eventualCollision(std::vector<glm::vec3> &pos,std::vector<glm::vec3> &vel){

			//Compute kd tree for the rigid body

			//treeNode* root=new  treeNode();
			//std::vector <glm::vec4> points;

			//for(int k=0;k<pos.size();k++)
			//{
			//	points.push_back(glm::vec4(pos[k],k));

			//}
			//root=SPHSimulation::computeKdTree(points,root,0);		

			//std::vector<treeNode*> broadlist;


			//SPHSimulation::broadPhase(root,broadlist,triList);


			//std::vector<glm::vec4>  ptl;
			//for(int i=0;i<broadlist.size();i++){

			//	glm::vec4 pt=broadlist[i]->location;

			//	for (int k = 0; k < broadlist[i]->points.size(); k++)
			//	{
			//		ptl.push_back(broadlist[i]->points[k]);
			//	}

			//}

			//for (int k = 0; k < ptl.size(); k++)
			//{
			//	for (int i = 0; i < triList.size(); i+=3)
			//	{
			//		glm::vec3 *T=new glm::vec3[3];;
			//		T=(glm::vec3*) malloc (3);
			//		T[0]=triList[i];
			//		T[1]=triList[i+1];
			//		T[2]=triList[i+2];
			//		glm::vec3 tp=glm::vec3(ptl[k].x,ptl[k].y,ptl[k].z);
			//		bool currentPos=SPHSimulation::sphereIntersectionTriangle3D(T,tp,globalRadius);
			//		glm::vec3 np=tp+vel[ptl[k].a]*dt;
			//		bool nextPos=SPHSimulation::sphereIntersectionTriangle3D(T,np,globalRadius);
			//		Segment R;
			//		R.P0=tp;
			//		R.P1=np;
			//		int lineInt=SPHSimulation::intersect3D_RayTriangle( R,T);

			//		if(currentPos||nextPos||lineInt){
			//			glm::vec3 V = glm::cross(T[1] - T[0], T[2] - T[0]);

			//			float angle=glm::angle(V,vel[ptl[k].a]);
			//			float ca=cos(angle);
			//			vel[ptl[k].a]=-ca*vel[ptl[k].a]+(1-ca)*vel[ptl[k].a];

			//		}

			//	}	
			//}
			for (int k = 0; k < pos.size(); k++)
			{
				for (int i = 0; i < triList.size(); i+=3)
				{
					glm::vec3 *T=new glm::vec3[3];

					T[0]=glm::vec3(triList[i].x,triList[i].y,triList[i].z);
					T[1]=glm::vec3(triList[i+1]);
					T[2]=glm::vec3(triList[i+2]);

					glm::vec3 tp=pos[k];

					bool currentPos=SPHSimulation::sphereIntersectionTriangle3D(T,tp,globalRadius);

					glm::vec3 np=tp+vel[k]*dt;
					bool nextPos=SPHSimulation::sphereIntersectionTriangle3D(T,np,globalRadius);

					Segment R;
					R.P0=tp;
					R.P1=np;

					bool lineInt=false; 
					lineInt=SPHSimulation::intersect3D_RayTriangle( R,T);

					if(currentPos||nextPos||lineInt){

						glm::vec3 V = glm::cross(T[1] - T[0], T[2] - T[0]);
						// cosθ = a.b / |a||b| 
						float an=glm::dot(V,vel[k])/(glm::length(V)*glm::length(vel[k]));
						an=std::abs(an);
						float sa=std::sqrt(1-(an*an));
						vel[k]=(-5.35f*an*vel[k]+(sa)*vel[k]);

						break;
					}

				}	
			}

		}
		//giving bullet phyiscs instructions what to do when a collision pair is detected
		// using a callback function on near phase detections

		//void myCdaCallback(btBroadphasePair& collisionPair,
		//	btCollisionDispatcher& dispatcher, btDispatcherInfo& dispatchInfo){


		//}

		////using the collision detection module of BUllet physics
		//void SPHSimulation::eventualCollision(std::vector<glm::vec3> &pos,std::vector<glm::vec3> &vel){


		//	//keep track of the shapes, we release memory at exit.
		//	//make sure to re-use collision shapes among rigid bodies whenever possible!

		//	btVector3 *positions;
		//	positions=(btVector3*)malloc(sizeof(float)*sphereCount*4);

		//	for (int i = 0; i < sphereCount; i++)
		//	{
		//		positions[i].setX((float)pos[i].x);
		//		positions[i].setY((float)pos[i].y);
		//		positions[i].setZ((float)pos[i].z);

		//	}

		//	btVector3 halfExt(1,1.51,1);
		//	btCollisionShape* box=new btBoxShape (halfExt);
		//	btCollisionShape* fallShape = new btSphereShape(globalRadius);



		//	btDefaultCollisionConfiguration* collisionConfiguration = new btDefaultCollisionConfiguration();
		//	btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collisionConfiguration);
		//	btVector3	worldAabbMin(-100,-100,-100);
		//	btVector3	worldAabbMax(100,100,100);

		//	btAxisSweep3*	broadphase = new btAxisSweep3(worldAabbMin,worldAabbMax);
		//	int numObj=SPHSimulation::sphereCount+1;
		//	btCollisionObject*	objects;

		//	btSequentialImpulseConstraintSolver* solver = new btSequentialImpulseConstraintSolver;
		//	btDiscreteDynamicsWorld* dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher,broadphase,solver,collisionConfiguration);
		//	dynamicsWorld->setGravity(btVector3(0,-10,0));

		//	btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0,1,0),0);

		//	btDefaultMotionState* groundMotionState = new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1),btVector3(0,0,-1)));
		//	btRigidBody::btRigidBodyConstructionInfo groundRigidBodyCI(0,groundMotionState,groundShape,btVector3(0,0,0));
		//	btRigidBody* groundRigidBody = new btRigidBody(groundRigidBodyCI);
		//	dynamicsWorld->addRigidBody(groundRigidBody);

		//	btDefaultMotionState* groundMotionState1 = new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1),btVector3(0,-1,0)));
		//	btRigidBody::btRigidBodyConstructionInfo groundRigidBodyCI1(0,groundMotionState1,groundShape,btVector3(0,0,0));
		//	btRigidBody* groundRigidBody1 = new btRigidBody(groundRigidBodyCI1);
		//	//dynamicsWorld->addRigidBody(groundRigidBody1);

		//	btDefaultMotionState* boxMotionState = new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1),btVector3(0,0,-1)));
		//	btRigidBody::btRigidBodyConstructionInfo boxRigidBodyCI(0,boxMotionState,box,btVector3(0,0,0));
		//	btRigidBody* boxRigidBody = new btRigidBody(boxRigidBodyCI);
		//	dynamicsWorld->addRigidBody(boxRigidBody);

		//	int bodiesBefore=2;

		//	std::vector<btRigidBody*> fallRigidBody;
		//	for (int i = 0; i < sphereCount; i++)
		//	{
		//		btRigidBody* frb;

		//		btDefaultMotionState* fallMotionState =
		//			new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1),positions[i]));
		//		btVector3 lv((float)vel[i].x,(float)vel[i].y,(float)vel[i].z);

		//		btRigidBody::btRigidBodyConstructionInfo fallRigidBodyCI(commonMass,fallMotionState,fallShape,fallInertia);
		//		frb=new btRigidBody(fallRigidBodyCI);
		//		frb->setLinearVelocity(lv);
		//		fallRigidBody.push_back(frb);
		//		dynamicsWorld->addRigidBody(frb);

		//	}

		//	dynamicsWorld->stepSimulation(dt,2);


		//	for (int i = 0; i < sphereCount; i++)
		//	{
		//		btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[i+bodiesBefore];
		//		btRigidBody* body = btRigidBody::upcast(obj);

		//		btTransform trans;
		//		fallRigidBody[i]->getMotionState()->getWorldTransform(trans);
		//		//btVector3 vtmp;
		//		btVector3 vtmp1;
		//		//vtmp=fallRigidBody[i-1]->getLinearVelocity();
		//		vtmp1=body->getLinearVelocity();
		//		pos[i]=glm::vec3( trans.getOrigin().getX(), trans.getOrigin().getY(), trans.getOrigin().getZ());
		//		vel[i]=glm::vec3(vtmp1.getX(),vtmp1.getY(),vtmp1.getZ());
		//		btVector3 aa;
		//	}

		//	delete fallShape;
		//	delete dynamicsWorld;
		//	delete solver;
		//	delete collisionConfiguration;
		//	delete dispatcher;
		//	delete broadphase;
		//	free(positions);

		//}
		void SPHSimulation::updateDT(float &dt,float vmax,float etaMax,float c, float h){
			//see equaqtion 13
			float t1=h/(vmax+c);
			float t2=(h*h)/(6*etaMax);
			if(t1<t2)
				dt=0.1*t1;
			else
				dt=0.1*t2;
		}
		void SPHSimulation::cudaTest(){
#ifdef ONGPU

			gpu.initializeGPU();
			std::cout<< std::endl<< std::endl<<"$$$$$$$$$$GPU INIT DONE$$$$$$$$$$$$$$$"<< std::endl<< std::endl;
#endif


		}

		bool SPHSimulation::loadOBJ(const char * path, 
			std::vector<glm::vec3> & out_vertices){

				printf("Loading OBJ file %s...\n", path);

				std::vector<unsigned int> vertexIndices;

				int count=0;
				FILE * file = fopen(path, "r");
				if( file == NULL ){
					printf("Impossible to open the file ! \n");
					return false;
				}

				while( 1 )
				{

					char lineHeader[128];
					// read the first word of the line
					int res = fscanf(file, "%s", lineHeader);
					if (res == EOF)
						break; // EOF = End Of File. Quit the loop.

					// else : parse lineHeader

					if ( strcmp( lineHeader, "v" ) == 0 ){
						glm::vec3 vertex;

						fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
						out_vertices.push_back(vertex);
						count++;
					}
					else{
						// Probably a comment, eat up the rest of the line
						char stupidBuffer[1000];
						fgets(stupidBuffer, 1000, file);
					}

					SPHSimulation::sphereCount=count;

				}


				return true;
		}

	}
}
