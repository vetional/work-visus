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
		const unsigned int SPHSimulation::sphereCount = 2000;
#else
		const unsigned int SPHSimulation::sphereCount = 200;
#endif

		/*
		* SPHSimulation::SPHSimulation
		*/
		SPHSimulation::SPHSimulation(void) : view::AnimDataModule(), getDataSlot("getData", "Gets the data from the data source") {
			SPHSimulation::initData(m,roh,proh,droh,D,d,eta,pos,vel,gradV,acc,pacc,S,gradS,P,gradP,gradW,pie,W);
#ifdef ONGPU
			cudaTest();
#endif			
			this->getDataSlot.SetCallback(moldyn::MultiParticleDataCall::ClassName(), moldyn::MultiParticleDataCall::FunctionName(0), &SPHSimulation::getDataCallback);
			this->getDataSlot.SetCallback(moldyn::MultiParticleDataCall::ClassName(), moldyn::MultiParticleDataCall::FunctionName(1), &SPHSimulation::getExtentCallback);
			this->MakeSlotAvailable(&this->getDataSlot);

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
		float h=3*SPHSimulation::globalRadius;
		float dt=0.001f;//time step
		float vmax=100;//for dt update
		float etaMax=3;//for dt update
		float alpha=0.5;//bulk viscosity
		float epsilon=0.1f;//for velocity correction
		float roh0=1;//initial density
		float c=300000;// in m/sec
		float n=0.5f;//refer equations in paper for eta
		float jumpN=0.5f;//jump number


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
		glm::vec3 z(0.0f,0.0f,0.0f);
		glm::mat3 zero(z,z,z);
		float ze=0.0f;
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


				glm::vec3 velocity(0.02,0.02,.5);
				for(int i = 0;i<SPHSimulation::sphereCount;i++){
					m.push_back(1);
					roh.push_back(1);
					proh.push_back(1);
					droh.push_back(1);
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
					pos.push_back(z);
				}
				float a1=-0.3;float b1=-.15;float c1=.5;
				int dc=0;
				int side=5;
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
							a1+=globalRadius*2.5f;
						}
						a1=-0.3;
						b1+=globalRadius*2.5f;

					}c1+=globalRadius*2.5f;b1=-.15;

				}


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
			SPHSimulation::correctV(vel,m,W,epsilon,roh,pos,h);//see equaqtion 11
#endif
#ifdef ONGPU

			gpu.theBigLoop(m,roh,proh,droh,D,d,eta,pos,vel,gradV,acc,pacc,S,gradS,P,
				gradP,gradW,pie,W, h,sphereCount,epsilon,c,roh0,alpha,jumpN,dt);
#endif
			SPHSimulation::eventualCollision();

			SPHSimulation::updateDT(dt,vmax,etaMax,c,h);//see equaqtion 13

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

				void SPHSimulation::findNeighbours(int index, std::vector<glm::vec3> &pos, 
			std::vector<int> &nxi,float h){

				for(int i =0; i <pos.size();i++){



					if(glm::distance(pos[index],pos[i])<h){

						if(index!=i){
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
					gradW[j].z= multiplier*(pos[i].y-pos[nxi[j]].z);
				}
		}


		void SPHSimulation::updategradV(std::vector<float> &m,std::vector<float> &roh,
			std::vector<glm::vec3> &vel,std::vector<glm::mat3> &gradV,
			std::vector<glm::vec3> &gradW, int index,std::vector<int> &nxi)
		{

			float gv=0;
		
			for(int j=0;j<nxi.size();j++)
			{
				
				glm::vec3 v=(m[nxi[j]]/roh[nxi[j]])*(vel[nxi[j]]-vel[index]);
				gradV[index]+=glm::outerProduct(v,gradW[j]);

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
				eta[index]=(1-exp)*((1/d11)+(1/d11));
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

					ptemp=(m[nxi[j]])*((P[nxi[j]]/std::pow(roh[nxi[j]],2))+P[index]/std::pow(roh[index],2));
					sum+=ptemp*gradW[j];

				}
				gradP[index]=sum;
		}

		void SPHSimulation::updateAcc(std::vector<glm::vec3>&gradP,std::vector<glm::vec3> &gradS,
			std::vector<glm::vec3>&pie,std::vector<glm::vec3> &acc,std::vector<glm::vec3> &pacc){

				glm::vec3 gr(0.0f,-10.0001f,0.0f);
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
			std::vector<float> &W,float epsilon,std::vector<float> &roh,
			std::vector<glm::vec3> &pos,float h){//see equaqtion 11	

				for(int i = 0;i<m.size();i++){

					std::vector<int> nxi;
					SPHSimulation::findNeighbours(i,pos,nxi,h);
					glm::vec3 sum(0.0f,0.0f,0.0f);
					for(int j=0;j<nxi.size();j++)
					{
						sum+=((2*m[nxi[j]])/(roh[i]+roh[nxi[j]]))*(vel[nxi[j]]-vel[i])*W[j];
					}		

					vel[i]+=epsilon*sum;
					pos[i]=pos[i]+vel[i]*dt;
				}
		}
		void SPHSimulation::eventualCollision(){}

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
#endif


		}
	}
}
