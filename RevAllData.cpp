//Compilation info
//g++ RevAllData.cpp -o DATA -lgsl -lgslcblas -lSDL -lSDL_gfx -lSDL_ttf -lSDL_image -lpng -lm -lgomp -fopenmp -std=c++11
//Movie maker for a folder with only all the pics:
//convert -delay 10 -loop 1 `ls -v` WT.gif

#include <omp.h>                //Paralelización
#include <vector>               //std::vector http://www.cplusplus.com/reference/std::vector/
#include <string>               //http://www.cplusplus.com/reference/cstring/
#include <cmath>                //Librería matemática de C http://www.cplusplus.com/reference/cmath/
#include <iostream>             //Imprimir en Pantalla http://www.cplusplus.com/reference/iostream/
#include <iomanip>				//Parametric manipulators http://www.cplusplus.com/reference/iomanip/
#include <fstream>		        //Manejo con archivos C++ http://www.cplusplus.com/reference/fstream/
#include <sstream>				//Stream class to operate on strings. http://www.cplusplus.com/reference/sstream/stringstream/
#include <gsl/gsl_rng.h>        //Librerías GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>
#include <stdlib.h>             //Usado para itoa http://www.cplusplus.com/reference/cstdlib/
#include <stdio.h>              //Input/Output operations in C http://www.cplusplus.com/reference/cstdio/
#include <GL/gl.h>
#include <SDL/SDL.h>            	//Librerías SDL para los dibujos
#include <SDL/SDL_gfxPrimitives.h>
#include <SDL/SDL_ttf.h>

using namespace std;
////////////////////Simulation parameters//////////////////////
int const Ttot=96;
double const AT=1.; 	  		    //Intervalo de tiempo en el que mostrar el filamento y las concentraciones
int const MED=150;     			   	//Número de realizaciones para obtener la media
int const NBch=15;     			   	//Número data batches para obtener la media
int const MaxThr=10;
int const INITIAL_N_CELLS=30; 	 	//Número inicial de células
int const SEMILLA=13;			    //Semilla del generador de números aleatorios

const int WINDOW_WIDTH=1900;
const int WINDOW_HEIGHT=1000;
int const HISTOGRAM_SIZE=27;
int const BORDER_SIZE=6;

double const At=0.001; 			    //0.0167;		//Time step size (h)

///////////////////////System parameters///////////////////////
//Crecimiento y división
double GROWTH_RATE=0.08;					//Growth rate  [1/h]
double const MAX_SIZE=4.;    		    	//Maximum cell size
const double NOISE_AMP_initial_size=0.1; 	//Amplitud del ruido de la diferencia de tamaño inicial y de las hijas
const double NOISE_AMP_max_size=0.1;     	//Amplitud del ruido del tamaño máximo
const double NOISE_AMP_growth_rate=0.1;  	//Amplitud del ruido de la velocidad de crecimiento
//Diferenciación
const double MAX_LEVEL=110.0;            	//Concentración para diferenciarse [nM]
const double NOISE_AMP_max_level=0.15;    	//Amplitud del ruido del nivel para convertirse en heterociste
const double T_MEAN=12.0;                	//Tiempo medio para decidir diferenciarse  [h]
const double T_MIN=5.0;                	//Tiempo minimo para decidir diferenciarse  [h]
//Protoheterocyst
double T_MEAN_PROTO;               			//Tiempo medio de maduración  [h]
//Ruidos dynamica
const double NOISE_AMP_initial_Con=0.15;   	//Amplitud del ruido de las concentraciones iniciales
const double NOISE_AMP=0.15;               	//Amplitud del ruido de la evolución de las concentraciones
//PatX 
const double PatX=0.1;
gsl_rng *r; 										//Global random generator

class CCell
{
	private:

	//////////////////////Equation coeficients/////////////////////
	//HetR
	double beta_r;		//HetR basal production rate [nM/h]
	double rho_r;     //HetR regulated production rate [nM/h]
	double alpha_r;   //HetR degradation rate [1/h]
	double alpha_d;   //HetR dimer mediated degradation rate [1/nM]
	double K_r;       //HetR-Promoter equilibrium constant [nM]
	double alpha_GFP;	//GFP degradation rate [1/h]
	//PatA
	double rho_a;     //PatA regulated production rate [nM/h]
	double alpha_a;   //PatA degradation rate [1/h]
	double t_a;			//PatA/HetF-HetR activation equilibrium constant [nM^2]
	//HetF
	double EC_f;    	//Efective interactive HetF concentration [1]
	//PatS
	double rho_s;     //PatS regulated production rate [nM/h]
	double alpha_s;   //PatS degradation rate [1/h]
	double c_s;       //PatS-Inhibitor conversion rate [1/h]
	//HetN
	double rho_n;     //HetN basal production rate in heterocysts [nM/h]
	double alpha_n;   //HetN degradation rate [1/h]
	double c_n; 	   //HetN-Inhibitor conversion rate [1/h]
	//Inhibitor
	double alpha_i;   //Inhibitor degradation rate [1/h]
	double d_i; 	   //Inhibitor difusion rate [1/h]
	double Kd;	      //Inhibitor-HetR equilibrium constant [nM]
	//Nitrogen
	double rho_NF;	   //Nitrogen basal production rate in heterocysts [nM/h]
	double alpha_NF;	//Nitrogen degradation rate [1/h]
	double d_NF;	   //Nitrogen difusion rate [1/h]
	double K_NF;		//Nitrogen-HetR equilibrium constant [nM]
	//Frontera
	double d_border;  //Modificador difusión en los bordes
	//Noise
	double CellNoise;//Noise in celular processes
	double RegNoise; //Noise in concentration regulatory processes
	///////////////////////////////////////////////////////////////

	bool vegetative, protoheterocist;
	double size, max_size, max_level, cumulative_HetR, PROTO_T, DifD_T;
	double HetR, PatA, PatS, HetN, Inhb, fixN, PatS_old, HetN_old, Inhb_old, fixN_old, GFP;

	double Rnoise_c_s, Rnoise_c_n, Rnoise_c_is, Rnoise_c_in, Rnoise_d_i, Rnoise_d_NF; //Ruidos acoplados entre celulas

	public:
	CCell(std::vector<double> Param)
	{
		beta_r=Param[0];
		rho_r=Param[1];
		alpha_r=Param[2];
		alpha_d=Param[3];
		K_r=Param[4];
		alpha_GFP=Param[5];
		rho_a=Param[6];
		alpha_a=Param[7];
		t_a=Param[8];
		EC_f=Param[9];
		rho_s=Param[10];
		alpha_s=Param[11];
		c_s=Param[12];
		rho_n=Param[13];
		alpha_n=Param[14];
		c_n=Param[15];
		alpha_i=Param[16];
		d_i=Param[17];
		Kd=Param[18];
		rho_NF=Param[19];
		alpha_NF=Param[20];
		d_NF=Param[21];
		K_NF=Param[22];
		d_border=Param[23];
		CellNoise=Param[24];
		RegNoise=Param[25];

		vegetative=true;
		protoheterocist=false;
		size=MAX_SIZE*0.75*abs(1.+CellNoise*NOISE_AMP_initial_size*gsl_ran_gaussian(r,1));
		HetR=beta_r/alpha_r*abs(1.+RegNoise*NOISE_AMP_initial_Con*gsl_ran_gaussian(r,1));
		PatA=0.;
		PatS=0.;
		HetN=0.;
		Inhb=0.;
		fixN=0.;
		Set_cell(0);			//Tamaño medio inicial 3 micras
	}

	double Get_size() {return size;}
	double GeT_MEAN_size() {return max_size;}
	double GeT_MEAN_level() {return max_level;}
	double Set_size_cell1()
	{
		size=size*0.5*(1.+CellNoise*NOISE_AMP_initial_size*gsl_ran_gaussian(r,1));
		return (max_size-size);
	}
	void Set_size_cell2(double size2) {size=size2;}
	void Set_cell(double OHetRA)
	{
		max_size=MAX_SIZE*abs(1.+CellNoise*NOISE_AMP_max_size*gsl_ran_gaussian(r,1)); //Establece el tamaño máximo intrínseco al que ha dividirse cada célula
		max_level=MAX_LEVEL*abs(1.+RegNoise*NOISE_AMP_max_level*gsl_ran_gaussian(r,1)); //Establece el nivel máximo de concentración acumulado para convertirse en heterociste
		cumulative_HetR=OHetRA;
		PROTO_T=0;
		DifD_T=0;
	}
	double Get_HetR() {return HetR;}
	double Get_GFP() {return GFP;}
	double Get_PatA() {return PatA;}
	double Get_PatS() {return PatS;}
	double Get_HetN() {return HetN;}
	double Get_Inhb() {return Inhb;}
	double Get_fixN() {return fixN;}
	double Get_cumulative_HetR() {return cumulative_HetR;}
	double Get_PROTO_T(){ return PROTO_T;}
	double Get_DifD_T(){ return DifD_T;}

	bool Get_vegetative() {return vegetative;}
	bool Get_protoheterocist() {return protoheterocist;}
	void Transform_protoheterocist() {protoheterocist=true;}
	void Transform_heterocist() {vegetative=false; protoheterocist=false;}
	void Advance_At_vegetative(std::vector<CCell> &filament, int position)
	{
		int n_cells=filament.size();  //Variable temporal guarda número de células del filamento
		double fhp,rhp;  //Regulation
		double C_PatS_L, C_PatS_R, C_IS_L, C_IS_R, C_IN_L, C_IN_R, D_Inhb_L, D_Inhb_R, D_fixN_L, D_fixN_R;  //Difusion
		double noise_t_F, noise_t_A, Lnoise_c_s, Lnoise_c_n, Lnoise_c_is, Lnoise_c_in, Lnoise_d_i, Lnoise_d_NF, noise_proR;		//Ruidos acoplados

		//Hay que hacer copias para actualizar correctamente la integral y para calcular bien el laplaciano
		PatS_old=PatS;
		HetN_old=HetN;
		Inhb_old=Inhb;
		fixN_old=fixN;

		if(position==0) //Difusion flujo en los extremos
		{
			//Conversión PatS-Inhb
			C_PatS_L=0;
			C_PatS_R=c_s*PatS_old;
			Lnoise_c_s=0;
			Rnoise_c_s=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=0;
			C_IS_R=c_s*(filament.at(position+1).PatS);
			Lnoise_c_is=0;
			Rnoise_c_is=gsl_ran_gaussian(r,1)*sqrt(C_IS_R/At);
			C_IN_L=0;
			C_IN_R=c_n*(filament.at(position+1).HetN);
			Lnoise_c_in=0;
			Rnoise_c_in=gsl_ran_gaussian(r,1)*sqrt(C_IN_R/At);
			//Difusion Inhb
			D_Inhb_L=d_i*d_border*Inhb_old;
			D_Inhb_R=d_i*(Inhb_old-filament.at(position+1).Inhb);
			Lnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(D_Inhb_L/At);
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(d_i*(Inhb_old+filament.at(position+1).Inhb)/At);
			//Difusion fixN
			D_fixN_L=d_NF*d_border*fixN_old;
			D_fixN_R=d_NF*(fixN_old-filament.at(position+1).fixN);
			Lnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(D_fixN_L/At);
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(d_NF*(fixN_old+filament.at(position+1).fixN)/At);
		}
		else if(position==(n_cells-1))
		{
			//Conversión PatS-Inhb
			C_PatS_L=c_s*PatS_old;
			C_PatS_R=0;
			Lnoise_c_s=filament.at(position-1).Rnoise_c_is;
			Rnoise_c_s=0;
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=c_s*(filament.at(position-1).PatS_old);
			C_IS_R=0;
			Lnoise_c_is=filament.at(position-1).Rnoise_c_s;
			Rnoise_c_is=0;
			C_IN_L=c_n*(filament.at(position-1).HetN_old);
			C_IN_R=0;
			Lnoise_c_in=filament.at(position-1).Rnoise_c_n;
			Rnoise_c_in=0;
			//Difusion Inhb
			D_Inhb_L=d_i*(Inhb_old-filament.at(position-1).Inhb_old);
			D_Inhb_R=d_i*d_border*Inhb_old;
			Lnoise_d_i=filament.at(position-1).Rnoise_d_i;
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(D_Inhb_R/At);
			//Difusion fixN
			D_fixN_L=d_NF*(fixN_old-filament.at(position-1).fixN_old);
			D_fixN_R=d_NF*d_border*fixN_old;
			Lnoise_d_NF=filament.at(position-1).Rnoise_d_NF;
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(D_fixN_R/At);
		}
		else
		{
			//Conversión PatS-Inhb
			C_PatS_L=c_s*PatS_old;
			C_PatS_R=c_s*PatS_old;
			Lnoise_c_s=filament.at(position-1).Rnoise_c_is;
			Rnoise_c_s=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=c_s*(filament.at(position-1).PatS_old);
			C_IS_R=c_s*(filament.at(position+1).PatS);
			Lnoise_c_is=filament.at(position-1).Rnoise_c_s;
			Rnoise_c_is=gsl_ran_gaussian(r,1)*sqrt(C_IS_R/At);
			C_IN_L=c_n*(filament.at(position-1).HetN_old);
			C_IN_R=c_n*(filament.at(position+1).HetN);
			Lnoise_c_in=filament.at(position-1).Rnoise_c_n;
			Rnoise_c_in=gsl_ran_gaussian(r,1)*sqrt(C_IN_R/At);
			//Difusion Inhb
			D_Inhb_L=d_i*(Inhb_old-filament.at(position-1).Inhb_old);
			D_Inhb_R=d_i*(Inhb_old-filament.at(position+1).Inhb);
			Lnoise_d_i=filament.at(position-1).Rnoise_d_i;
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(d_i*(Inhb_old+filament.at(position+1).Inhb)/At);
			//Difusion fixN
			D_fixN_L=d_NF*(fixN_old-filament.at(position-1).fixN_old);
			D_fixN_R=d_NF*(fixN_old-filament.at(position+1).fixN);
			Lnoise_d_NF=filament.at(position-1).Rnoise_d_NF;
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(d_NF*(fixN_old+filament.at(position+1).fixN)/At);
		}

		fhp=(gsl_pow_2(HetR/K_r)*((EC_f*PatA/t_a)+EC_f))/(1.+(gsl_pow_2(HetR/K_r)*((EC_f*PatA/t_a)+EC_f))+gsl_pow_2(Inhb_old/Kd)+fixN_old/K_NF);
		//~ rhp=(gsl_pow_4(1.9*HetR/K_r)*((EC_f*PatA/t_a)+EC_f))/(1.+(gsl_pow_4(1.9*HetR/K_r)*((EC_f*PatA/t_a)+EC_f))+gsl_pow_4(Inhb_old/Kd)+fixN_old/K_NF);
		//~ noise_proR=RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((beta_r+rho_r*rhp)/At);
		//~ HetR+=At*(beta_r+rho_r*rhp-alpha_r*HetR*(1.+2.*alpha_d*HetR)+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((alpha_r*HetR*(1.+4.*alpha_d*HetR))/At)+noise_proR);
		
		noise_proR=RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((beta_r+rho_r*fhp)/At);
		HetR+=At*(beta_r+rho_r*fhp-alpha_r*HetR*(1.+2.*alpha_d*HetR)+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((alpha_r*HetR*(1.+4.*alpha_d*HetR))/At)+noise_proR);
		if(HetR<0) HetR=0;

		//~ GFP+=At*(beta_r+rho_r*fhp-alpha_GFP*GFP+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((alpha_GFP*GFP)/At)+noise_proR);
		//~ if(GFP<0) GFP=0;

		PatA+=At*(rho_a*fhp-alpha_a*PatA+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((rho_a*fhp+alpha_a*PatA)/At));
		if(PatA<0) PatA=0;

		PatS+=At*(rho_s*fhp-C_PatS_L-C_PatS_R-alpha_s*PatS_old+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((rho_s*fhp+alpha_s*PatS_old)/At)+Lnoise_c_s+Rnoise_c_s));
		if(PatS<0) PatS=0;

		Inhb+=At*(C_IS_L+C_IS_R+C_IN_L+C_IN_R-D_Inhb_L-D_Inhb_R-alpha_i*Inhb_old+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((alpha_i*Inhb_old)/At)+Lnoise_c_is+Rnoise_c_is+Lnoise_c_in+Rnoise_c_in+Lnoise_d_i+Rnoise_d_i));
		if(Inhb<0) Inhb=0;

		fixN+=At*(-alpha_NF*fixN_old-D_fixN_L-D_fixN_R+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((alpha_NF*fixN_old)/At)+Lnoise_d_NF+Rnoise_d_NF));
		if(fixN<0) fixN=0;

		if(HetR>max_level)
		{
			cumulative_HetR+=HetR*At;
			DifD_T+=At;
		}
		else
		{
			cumulative_HetR=0;
			DifD_T=0;
		}

		size+=At*(GROWTH_RATE+CellNoise*NOISE_AMP_growth_rate*gsl_ran_gaussian(r,1)*sqrt(GROWTH_RATE/At)); 		//Crece la célula
	}

	void Advance_At_protoheterocist(std::vector<CCell> &filament, int position)
	{
		int n_cells=filament.size();  //Variable temporal guarda número de células del filamento
		double fhp;  //Regulation
		double C_PatS_L, C_PatS_R, C_HetN_L, C_HetN_R, C_IS_L, C_IS_R, C_IN_L, C_IN_R, D_Inhb_L, D_Inhb_R, D_fixN_L, D_fixN_R;  //Difusion
		double noise_t_F, noise_t_A, Lnoise_c_s, Lnoise_c_n, Lnoise_c_is, Lnoise_c_in, Lnoise_d_i, Lnoise_d_NF, noise_proR;		//Ruidos acoplados

		//Hay que hacer copias para actualizar correctamente la integral y para calcular bien el laplaciano
		PatS_old=PatS;
		HetN_old=HetN;
		Inhb_old=Inhb;
		fixN_old=fixN;

		if(position==0) //Difusion flujo en los extremos
		{
			//Conversión PatS-Inhb
			C_PatS_L=0;
			C_PatS_R=c_s*PatS_old;
			Lnoise_c_s=0;
			Rnoise_c_s=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			//Conversión HetN-Inhb
			C_HetN_L=0;
			C_HetN_R=c_n*HetN_old;
			Lnoise_c_n=0;
			Rnoise_c_n=gsl_ran_gaussian(r,1)*sqrt(C_HetN_R/At);
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=0;
			C_IS_R=c_s*(filament.at(position+1).PatS);
			Lnoise_c_is=0;
			Rnoise_c_is=gsl_ran_gaussian(r,1)*sqrt(C_IS_R/At);
			C_IN_L=0;
			C_IN_R=c_n*(filament.at(position+1).HetN);
			Lnoise_c_in=0;
			Rnoise_c_in=gsl_ran_gaussian(r,1)*sqrt(C_IN_R/At);
			//Difusion Inhb
			D_Inhb_L=d_i*d_border*Inhb_old;
			D_Inhb_R=d_i*(Inhb_old-filament.at(position+1).Inhb);
			Lnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(D_Inhb_L/At);
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(d_i*(Inhb_old+filament.at(position+1).Inhb)/At);
			//Difusion fixN
			D_fixN_L=d_NF*d_border*fixN_old;
			D_fixN_R=d_NF*(fixN_old-filament.at(position+1).fixN);
			Lnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(D_fixN_L/At);
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(d_NF*(fixN_old+filament.at(position+1).fixN)/At);
		}
		else if(position==(n_cells-1))
		{
			//Conversión PatS-Inhb
			C_PatS_L=c_s*PatS_old;
			C_PatS_R=0;
			Lnoise_c_s=filament.at(position-1).Rnoise_c_is;
			Rnoise_c_s=0;
			//Conversión HetN-Inhb
			C_HetN_L=c_n*HetN_old;
			C_HetN_R=0;
			Lnoise_c_n=filament.at(position-1).Rnoise_c_in;
			Rnoise_c_n=0;
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=c_s*(filament.at(position-1).PatS_old);
			C_IS_R=0;
			Lnoise_c_is=filament.at(position-1).Rnoise_c_s;
			Rnoise_c_is=0;
			C_IN_L=c_n*(filament.at(position-1).HetN_old);
			C_IN_R=0;
			Lnoise_c_in=filament.at(position-1).Rnoise_c_n;
			Rnoise_c_in=0;
			//Difusion Inhb
			D_Inhb_L=d_i*(Inhb_old-filament.at(position-1).Inhb_old);
			D_Inhb_R=d_i*d_border*Inhb_old;
			Lnoise_d_i=filament.at(position-1).Rnoise_d_i;
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(D_Inhb_R/At);
			//Difusion fixN
			D_fixN_L=d_NF*(fixN_old-filament.at(position-1).fixN_old);
			D_fixN_R=d_NF*d_border*fixN_old;
			Lnoise_d_NF=filament.at(position-1).Rnoise_d_NF;
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(D_fixN_R/At);
		}
		else
		{
			//Conversión PatS-Inhb
			C_PatS_L=c_s*PatS_old;
			C_PatS_R=c_s*PatS_old;
			Lnoise_c_s=filament.at(position-1).Rnoise_c_is;
			Rnoise_c_s=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			//Conversión HetN-Inhb
			C_HetN_L=c_n*HetN_old;
			C_HetN_R=c_n*HetN_old;
			Lnoise_c_n=filament.at(position-1).Rnoise_c_in;
			Rnoise_c_n=gsl_ran_gaussian(r,1)*sqrt(C_HetN_R/At);
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=c_s*(filament.at(position-1).PatS_old);
			C_IS_R=c_s*(filament.at(position+1).PatS);
			Lnoise_c_is=filament.at(position-1).Rnoise_c_s;
			Rnoise_c_is=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			C_IN_L=c_n*(filament.at(position-1).HetN_old);
			C_IN_R=c_n*(filament.at(position+1).HetN);
			Lnoise_c_in=filament.at(position-1).Rnoise_c_n;
			Rnoise_c_in=gsl_ran_gaussian(r,1)*sqrt(C_HetN_R/At);
			//Difusion Inhb
			D_Inhb_L=d_i*(Inhb_old-filament.at(position-1).Inhb_old);
			D_Inhb_R=d_i*(Inhb_old-filament.at(position+1).Inhb);
			Lnoise_d_i=filament.at(position-1).Rnoise_d_i;
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(d_i*(Inhb_old+filament.at(position+1).Inhb)/At);
			//Difusion fixN
			D_fixN_L=d_NF*(fixN_old-filament.at(position-1).fixN_old);
			D_fixN_R=d_NF*(fixN_old-filament.at(position+1).fixN);
			Lnoise_d_NF=filament.at(position-1).Rnoise_d_NF;
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(d_NF*(fixN_old+filament.at(position+1).fixN)/At);
		}

		fhp=(gsl_pow_2(HetR/K_r)*((EC_f*PatA/t_a)+EC_f))/(1.+(gsl_pow_2(HetR/K_r)*((EC_f*PatA/t_a)+EC_f))+gsl_pow_2(Inhb_old/Kd)+fixN_old/K_NF);

		//~ noise_proR=RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((beta_r+rho_r*fhp)/At);
		//~ HetR+=At*(beta_r+rho_r*fhp-alpha_r*HetR*(1.+2.*alpha_d*HetR)+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((alpha_r*HetR*(1.+4.*alpha_d*HetR))/At)+noise_proR);
		//~ if(HetR<0) HetR=0;

		//~ GFP+=At*(beta_r+rho_r*fhp-alpha_GFP*GFP+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((alpha_GFP*GFP)/At)+noise_proR);
		//~ if(GFP<0) GFP=0;

		//~ PatA+=At*(rho_a*fhp-alpha_a*PatA+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((rho_a*fhp+alpha_a*PatA)/At));
		//~ if(PatA<0) PatA=0;

		PatS+=At*(rho_s*fhp*(1-(PROTO_T/T_MEAN_PROTO))-C_PatS_L-C_PatS_R-alpha_s*PatS_old+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((rho_s*fhp*(1-(PROTO_T/T_MEAN_PROTO))+alpha_s*PatS_old)/At)+Lnoise_c_s+Rnoise_c_s));
		if(PatS<0) PatS=0;

		HetN+=At*(rho_n*(PROTO_T/T_MEAN_PROTO)-alpha_n*HetN_old-C_HetN_L-C_HetN_R+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((rho_n*(PROTO_T/T_MEAN_PROTO)+alpha_n*HetN_old)/At)+Lnoise_c_n+Rnoise_c_n));
		if(HetN<0) HetN=0;

		Inhb+=At*(C_IS_L+C_IS_R+C_IN_L+C_IN_R-D_Inhb_L-D_Inhb_R-alpha_i*Inhb_old+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((alpha_i*Inhb_old)/At)+Lnoise_c_is+Rnoise_c_is+Lnoise_c_in+Rnoise_c_in+Lnoise_d_i+Rnoise_d_i));
		if(Inhb<0) Inhb=0;

		fixN+=At*(rho_NF-alpha_NF*fixN_old-D_fixN_L-D_fixN_R+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((rho_NF+alpha_NF*fixN_old)/At)+Lnoise_d_NF+Rnoise_d_NF));
		if(fixN<0) fixN=0;

		PROTO_T+=At;

		if(size<max_size) size+=At*(GROWTH_RATE+CellNoise*NOISE_AMP_growth_rate*gsl_ran_gaussian(r,1)*sqrt(GROWTH_RATE/At)); 		//Crece la célula heterociste
	}

	void Advance_At_heterocist(std::vector<CCell> &filament, int position)
	{
		int n_cells=filament.size();  //Variable temporal guarda número de células del filamento
		double fhp;  //Regulation
		double C_PatS_L, C_PatS_R, C_HetN_L, C_HetN_R, C_IS_L, C_IS_R, C_IN_L, C_IN_R, D_Inhb_L, D_Inhb_R, D_fixN_L, D_fixN_R;  //Difusion
		double noise_t_F, noise_t_A, Lnoise_c_s, Lnoise_c_n, Lnoise_c_is, Lnoise_c_in, Lnoise_d_i, Lnoise_d_NF, noise_proR;		//Ruidos acoplados

		//Hay que hacer copias para actualizar correctamente la integral y para calcular bien el laplaciano
		PatS_old=PatS;
		HetN_old=HetN;
		Inhb_old=Inhb;
		fixN_old=fixN;

		if(position==0) //Difusion flujo en los extremos
		{
			//Conversión PatS-Inhb
			C_PatS_L=0;
			C_PatS_R=c_s*PatS_old;
			Lnoise_c_s=0;
			Rnoise_c_s=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			//Conversión HetN-Inhb
			C_HetN_L=0;
			C_HetN_R=c_n*HetN_old;
			Lnoise_c_n=0;
			Rnoise_c_n=gsl_ran_gaussian(r,1)*sqrt(C_HetN_R/At);
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=0;
			C_IS_R=c_s*(filament.at(position+1).PatS);
			Lnoise_c_is=0;
			Rnoise_c_is=gsl_ran_gaussian(r,1)*sqrt(C_IS_R/At);
			C_IN_L=0;
			C_IN_R=c_n*(filament.at(position+1).HetN);
			Lnoise_c_in=0;
			Rnoise_c_in=gsl_ran_gaussian(r,1)*sqrt(C_IN_R/At);
			//Difusion Inhb
			D_Inhb_L=d_i*d_border*Inhb_old;
			D_Inhb_R=d_i*(Inhb_old-filament.at(position+1).Inhb);
			Lnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(D_Inhb_L/At);
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(d_i*(Inhb_old+filament.at(position+1).Inhb)/At);
			//Difusion fixN
			D_fixN_L=d_NF*d_border*fixN_old;
			D_fixN_R=d_NF*(fixN_old-filament.at(position+1).fixN);
			Lnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(D_fixN_L/At);
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(d_NF*(fixN_old+filament.at(position+1).fixN)/At);
		}
		else if(position==(n_cells-1))
		{
			//Conversión PatS-Inhb
			C_PatS_L=c_s*PatS_old;
			C_PatS_R=0;
			Lnoise_c_s=filament.at(position-1).Rnoise_c_is;
			Rnoise_c_s=0;
			//Conversión HetN-Inhb
			C_HetN_L=c_n*HetN_old;
			C_HetN_R=0;
			Lnoise_c_n=filament.at(position-1).Rnoise_c_in;
			Rnoise_c_n=0;
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=c_s*(filament.at(position-1).PatS_old);
			C_IS_R=0;
			Lnoise_c_is=filament.at(position-1).Rnoise_c_s;
			Rnoise_c_is=0;
			C_IN_L=c_n*(filament.at(position-1).HetN_old);
			C_IN_R=0;
			Lnoise_c_in=filament.at(position-1).Rnoise_c_n;
			Rnoise_c_in=0;
			//Difusion Inhb
			D_Inhb_L=d_i*(Inhb_old-filament.at(position-1).Inhb_old);
			D_Inhb_R=d_i*d_border*Inhb_old;
			Lnoise_d_i=filament.at(position-1).Rnoise_d_i;
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(D_Inhb_R/At);
			//Difusion fixN
			D_fixN_L=d_NF*(fixN_old-filament.at(position-1).fixN_old);
			D_fixN_R=d_NF*d_border*fixN_old;
			Lnoise_d_NF=filament.at(position-1).Rnoise_d_NF;
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(D_fixN_R/At);
		}
		else
		{
			//Conversión PatS-Inhb
			C_PatS_L=c_s*PatS_old;
			C_PatS_R=c_s*PatS_old;
			Lnoise_c_s=filament.at(position-1).Rnoise_c_is;
			Rnoise_c_s=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			//Conversión HetN-Inhb
			C_HetN_L=c_n*HetN_old;
			C_HetN_R=c_n*HetN_old;
			Lnoise_c_n=filament.at(position-1).Rnoise_c_in;
			Rnoise_c_n=gsl_ran_gaussian(r,1)*sqrt(C_HetN_R/At);
			//Entrada del Inhb creado por las celulas adyacentes
			C_IS_L=c_s*(filament.at(position-1).PatS_old);
			C_IS_R=c_s*(filament.at(position+1).PatS);
			Lnoise_c_is=filament.at(position-1).Rnoise_c_s;
			Rnoise_c_is=gsl_ran_gaussian(r,1)*sqrt(C_PatS_R/At);
			C_IN_L=c_n*(filament.at(position-1).HetN_old);
			C_IN_R=c_n*(filament.at(position+1).HetN);
			Lnoise_c_in=filament.at(position-1).Rnoise_c_n;
			Rnoise_c_in=gsl_ran_gaussian(r,1)*sqrt(C_HetN_R/At);
			//Difusion Inhb
			D_Inhb_L=d_i*(Inhb_old-filament.at(position-1).Inhb_old);
			D_Inhb_R=d_i*(Inhb_old-filament.at(position+1).Inhb);
			Lnoise_d_i=filament.at(position-1).Rnoise_d_i;
			Rnoise_d_i=gsl_ran_gaussian(r,1)*sqrt(d_i*(Inhb_old+filament.at(position+1).Inhb)/At);
			//Difusion fixN
			D_fixN_L=d_NF*(fixN_old-filament.at(position-1).fixN_old);
			D_fixN_R=d_NF*(fixN_old-filament.at(position+1).fixN);
			Lnoise_d_NF=filament.at(position-1).Rnoise_d_NF;
			Rnoise_d_NF=gsl_ran_gaussian(r,1)*sqrt(d_NF*(fixN_old+filament.at(position+1).fixN)/At);
		}

		//~ fhp=(gsl_pow_2(HetR/K_r)*((EC_f*PatA/t_a)+EC_f))/(1.+(gsl_pow_2(HetR/K_r)*((EC_f*PatA/t_a)+EC_f))+gsl_pow_2(Inhb_old/Kd)+fixN_old/K_NF);

		//~ noise_proR=RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((beta_r+rho_r*fhp)/At);
		//~ HetR+=At*(beta_r+rho_r*fhp-alpha_r*HetR*(1.+2.*alpha_d*HetR)+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((alpha_r*HetR*(1.+4.*alpha_d*HetR))/At)+noise_proR);
		//~ if(HetR<0) HetR=0;

		//~ GFP+=At*(beta_r+rho_r*fhp-alpha_GFP*GFP+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((alpha_GFP*GFP)/At)+noise_proR);
		//~ if(GFP<0) GFP=0;

		//~ PatA+=At*(rho_a*fhp-alpha_a*PatA+RegNoise*NOISE_AMP*gsl_ran_gaussian(r,1)*sqrt((rho_a*fhp+alpha_a*PatA)/At));
		//~ if(PatA<0) PatA=0;

		PatS+=At*(-C_PatS_L-C_PatS_R-alpha_s*PatS_old+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((alpha_s*PatS_old)/At)+Lnoise_c_s+Rnoise_c_s));
		if(PatS<0) PatS=0;

		HetN+=At*(rho_n-alpha_n*HetN_old-C_HetN_L-C_HetN_R+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((rho_n+alpha_n*HetN_old)/At)+Lnoise_c_n+Rnoise_c_n));
		if(HetN<0) HetN=0;

		Inhb+=At*(C_IS_L+C_IS_R+C_IN_L+C_IN_R-D_Inhb_L-D_Inhb_R-alpha_i*Inhb_old+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((alpha_i*Inhb_old)/At)+Lnoise_c_is+Rnoise_c_is+Lnoise_c_in+Rnoise_c_in+Lnoise_d_i+Rnoise_d_i));
		if(Inhb<0) Inhb=0;

		fixN+=At*(rho_NF-alpha_NF*fixN_old-D_fixN_L-D_fixN_R+RegNoise*NOISE_AMP*(gsl_ran_gaussian(r,1)*sqrt((rho_NF+alpha_NF*fixN_old)/At)+Lnoise_d_NF+Rnoise_d_NF));
		if(fixN<0) fixN=0;

		if(size<max_size) size+=At*(GROWTH_RATE+CellNoise*NOISE_AMP_growth_rate*gsl_ran_gaussian(r,1)*sqrt(GROWTH_RATE/At)); 		//Crece la célula heterociste
	}
};

std::string DoubleToString(double value, int precision)
{
	std::ostringstream streamValue;
	streamValue << std::setprecision(precision);
	streamValue << value;
	std::string strValue = streamValue.str();
	return strValue;
}

class CHistogram
{
	private:
	int array_size;
	int bord_array;
	int batch_array_size;
	int batch_bord_array;	
	
	double *histogram_veg_interval_med;
	double *histogram_veg_interval_disp;
	double *n_het_bound_hist_med;
	double *n_het_bound_hist_disp;
	double *n_cells_med;
	double *n_cells_disp;
	double *n_het_med;
	double *n_het_disp;
	double *n_het_boundaries_med;
	double *n_het_boundaries_disp;
	double *n_veg_med;
	double *n_veg_disp;
	
	double *clushet_size_med;
	double *clushet_size_disp;
	double *n_clushet_med;
	double *n_clushet_disp;
	
	double *p_het_med;
	double *p_het_disp;
	double *p_clushet_med;
	double *p_clushet_disp;
	
	
	double *concentration_HetR_filament_med;	
	//~ double **GFP_ConPos_n;
	//~ double **GFP_ConPos_med;
	//~ double **GFP_ConPos_disp;
	//~ double **HetR_ConPos_n;
	//~ double **HetR_ConPos_med;
	//~ double **HetR_ConPos_disp;

	//~ double **TGFP_ConPos_n;
	//~ double **TGFP_ConPos_med;
	//~ double **TGFP_ConPos_disp;
	//~ double **THetR_ConPos_n;
	//~ double **THetR_ConPos_med;
	//~ double **THetR_ConPos_disp;
	
	double *histogram_veg_interval_batch_med;
	double *n_het_bound_hist_batch_med;
	double *n_cells_batch_med;
	double *n_het_batch_med;
	double *n_het_boundaries_batch_med;
	
	double *histogram_veg_interval_batch_disp;
	double *n_het_bound_hist_batch_disp;
	double *n_cells_batch_disp;
	double *n_het_batch_disp;
	double *n_het_boundaries_batch_disp;
		
	double *clushet_size_batch_med;
	double *clushet_size_batch_disp;
	double *n_clushet_batch_med;
	double *n_clushet_batch_disp;
		
	double *p_het_batch_med;
	double *p_het_batch_disp;
	double *p_clushet_batch_med;
	double *p_clushet_batch_disp;

	public:
	CHistogram()
	{
		array_size=(int)HISTOGRAM_SIZE*((int)Ttot/AT+1);  //calcula el tamaño total del array dependiendo del tamaño del histograma y el espaciado de tiempo y tiempo total de simulación
		bord_array=BORDER_SIZE*((int)Ttot/AT+1);
		
		histogram_veg_interval_med=new double [array_size];
		histogram_veg_interval_disp=new double [array_size];		
		n_het_bound_hist_med=new double [bord_array];
		n_het_bound_hist_disp=new double [bord_array];	
		
		clushet_size_med=new double [bord_array];
		clushet_size_disp=new double [bord_array];
		
		histogram_veg_interval_batch_disp=new double [array_size];
		n_het_bound_hist_batch_disp=new double [bord_array];
		
		clushet_size_batch_disp=new double [bord_array];
		
		batch_array_size=NBch*(int)HISTOGRAM_SIZE*((int)Ttot/AT+1);
		batch_bord_array=NBch*BORDER_SIZE*((int)Ttot/AT+1);
		
		histogram_veg_interval_batch_med=new double [batch_array_size];
		n_het_bound_hist_batch_med=new double [batch_bord_array];
		
		clushet_size_batch_med=new double [batch_bord_array];
		

		for(long int i=0;i<batch_array_size;i++)
		{			
			histogram_veg_interval_batch_med[i]=0;	
				
			if(i<array_size)
			{
				histogram_veg_interval_med[i]=0;
				histogram_veg_interval_disp[i]=0;
				histogram_veg_interval_batch_disp[i]=0;	
			
				if(i<bord_array)
				{
					n_het_bound_hist_med[i]=0;
					n_het_bound_hist_disp[i]=0;
					n_het_bound_hist_batch_disp[i]=0;
					
					clushet_size_med[i]=0;
					clushet_size_disp[i]=0;
					clushet_size_batch_disp[i]=0;
				}
			}
						
			if(i<batch_bord_array)
			{
				n_het_bound_hist_batch_med[i]=0;
				clushet_size_batch_med[i]=0;				
			}		
		}

		n_cells_med=new double [array_size/HISTOGRAM_SIZE];
		n_cells_disp=new double [array_size/HISTOGRAM_SIZE];
		n_het_med=new double [array_size/HISTOGRAM_SIZE];
		n_het_disp=new double [array_size/HISTOGRAM_SIZE];
		n_veg_med=new double [array_size/HISTOGRAM_SIZE];
		n_veg_disp=new double [array_size/HISTOGRAM_SIZE];
		n_het_boundaries_med=new double [array_size/HISTOGRAM_SIZE];
		n_het_boundaries_disp=new double [array_size/HISTOGRAM_SIZE];
					
		n_clushet_med=new double [bord_array/BORDER_SIZE];
		n_clushet_disp=new double [bord_array/BORDER_SIZE];
		
		p_het_med=new double [array_size/HISTOGRAM_SIZE];
		p_het_disp=new double [array_size/HISTOGRAM_SIZE];	
		p_clushet_med=new double [bord_array/BORDER_SIZE];
		p_clushet_disp=new double [bord_array/BORDER_SIZE];
		
		concentration_HetR_filament_med=new double [array_size/HISTOGRAM_SIZE];
		
		n_cells_batch_disp=new double [array_size/HISTOGRAM_SIZE];
		n_het_batch_disp=new double [array_size/HISTOGRAM_SIZE];
		n_het_boundaries_batch_disp=new double [array_size/HISTOGRAM_SIZE];
				
		n_cells_batch_med=new double [batch_array_size/HISTOGRAM_SIZE];
		n_het_batch_med=new double [batch_array_size/HISTOGRAM_SIZE];
		n_het_boundaries_batch_med=new double [batch_array_size/HISTOGRAM_SIZE];
		
		n_clushet_batch_disp=new double [bord_array/BORDER_SIZE];
		n_clushet_batch_med=new double [batch_bord_array/BORDER_SIZE];
		
		p_het_batch_disp=new double [array_size/HISTOGRAM_SIZE];
		p_het_batch_med=new double [batch_array_size/HISTOGRAM_SIZE];
		p_clushet_batch_disp=new double [bord_array/BORDER_SIZE];
		p_clushet_batch_med=new double [batch_bord_array/BORDER_SIZE];

		
		for(long int i=0;i<batch_array_size/HISTOGRAM_SIZE;i++)
		{
			n_cells_batch_med[i]=0;
			n_het_batch_med[i]=0;
			n_het_boundaries_batch_med[i]=0;
			
			p_het_batch_med[i]=0;
			
			if(i<array_size/HISTOGRAM_SIZE)
			{
				n_cells_med[i]=0;
				n_cells_disp[i]=0;
				n_het_med[i]=0;
				n_het_disp[i]=0;
				n_veg_med[i]=0;
				n_veg_disp[i]=0;
				n_het_boundaries_med[i]=0;
				n_het_boundaries_disp[i]=0;
				concentration_HetR_filament_med[i]=0;
				
				p_het_med[i]=0;
				p_het_disp[i]=0;
				
				n_cells_batch_disp[i]=0;
				n_het_batch_disp[i]=0;
				n_het_boundaries_batch_disp[i]=0;
				
				p_het_batch_disp[i]=0;
			}
						
			if(i<bord_array/BORDER_SIZE)
			{
				n_clushet_med[i]=0;
				n_clushet_disp[i]=0;
				n_clushet_batch_disp[i]=0;	
				
				p_clushet_med[i]=0;
				p_clushet_disp[i]=0;
				p_clushet_batch_disp[i]=0;					
			}
			
			if(i<batch_bord_array/BORDER_SIZE)
			{
				n_clushet_batch_med[i]=0;
				
				p_clushet_batch_med[i]=0;
			}
		}

		//~ GFP_ConPos_n=new double* [(int)(Ttot/AT+1)];
		//~ GFP_ConPos_med=new double* [(int)(Ttot/AT+1)];
		//~ GFP_ConPos_disp=new double* [(int)(Ttot/AT+1)];
		//~ HetR_ConPos_n=new double* [(int)(Ttot/AT+1)];
		//~ HetR_ConPos_med=new double* [(int)(Ttot/AT+1)];
		//~ HetR_ConPos_disp=new double* [(int)(Ttot/AT+1)];

		//~ TGFP_ConPos_n=new double* [(int)(Ttot/AT+1)];
		//~ TGFP_ConPos_med=new double* [(int)(Ttot/AT+1)];
		//~ TGFP_ConPos_disp=new double* [(int)(Ttot/AT+1)];
		//~ THetR_ConPos_n=new double* [(int)(Ttot/AT+1)];
		//~ THetR_ConPos_med=new double* [(int)(Ttot/AT+1)];
		//~ THetR_ConPos_disp=new double* [(int)(Ttot/AT+1)];

		//~ for(long int i=0;i<(int)(Ttot/AT+1);i++)
		//~ {
			//~ GFP_ConPos_n[i]=new double [PMAX+1];
			//~ GFP_ConPos_med[i]=new double [PMAX+1];
			//~ GFP_ConPos_disp[i]=new double [PMAX+1];
			//~ HetR_ConPos_n[i]=new double [PMAX+1];
			//~ HetR_ConPos_med[i]=new double [PMAX+1];
			//~ HetR_ConPos_disp[i]=new double [PMAX+1];

			//~ TGFP_ConPos_n[i]=new double [TPMAX+1];
			//~ TGFP_ConPos_med[i]=new double [TPMAX+1];
			//~ TGFP_ConPos_disp[i]=new double [TPMAX+1];
			//~ THetR_ConPos_n[i]=new double [TPMAX+1];
			//~ THetR_ConPos_med[i]=new double [TPMAX+1];
			//~ THetR_ConPos_disp[i]=new double [TPMAX+1];
		//~ }

		//~ for(long int i=0;i<(int)(Ttot/AT+1);i++)
		//~ {
			//~ for(long int j=0;j<PMAX+1;j++)
			//~ {
				//~ GFP_ConPos_n[i][j]=0;
				//~ GFP_ConPos_med[i][j]=0;
				//~ GFP_ConPos_disp[i][j]=0;
				//~ HetR_ConPos_n[i][j]=0;
				//~ HetR_ConPos_med[i][j]=0;
				//~ HetR_ConPos_disp[i][j]=0;
			//~ }
			//~ for(long int j=0;j<TPMAX+1;j++)
			//~ {
				//~ TGFP_ConPos_n[i][j]=0;
				//~ TGFP_ConPos_med[i][j]=0;
				//~ TGFP_ConPos_disp[i][j]=0;
				//~ THetR_ConPos_n[i][j]=0;
				//~ THetR_ConPos_med[i][j]=0;
				//~ THetR_ConPos_disp[i][j]=0;
			//~ }
		//~ }
	}

	void Sum_histogram(std::vector<CCell> &filament, int indice_nt, int kSim)
	{
		int pos_het=0, temp_ncells=0, temp_nhet=0, temp_nbhet=0, fullH=0, clus_nhet=0, temp_nchet=0;
		int temp_histogram[HISTOGRAM_SIZE]={0}, temp_Bhist[BORDER_SIZE]={0}, temp_chet[BORDER_SIZE]={0};
		double temp_concentration_HetR=0;
				
		for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)   //El "for" lo avanza antes. Está bien
		{
			if(!(*it).Get_vegetative())
			{
				if(temp_nhet>=1 && (it-filament.begin()-pos_het-1)<HISTOGRAM_SIZE-1)  temp_histogram[(it-filament.begin()-pos_het-1)]++; //Si la distancia es mas pequeña que el tamaño-1 añade
				if(temp_nhet>=1 && (it-filament.begin()-pos_het-1)>=HISTOGRAM_SIZE-1) temp_histogram[HISTOGRAM_SIZE-1]++; //Si la distancia es mas grande que el tamaño-1 añade a la última posición del array
				pos_het=it-filament.begin();
				temp_nhet++;
				if(clus_nhet==0) temp_nchet++;	
				clus_nhet++;
			}
			else
			{
				if (clus_nhet!=0)
				{
					temp_chet[clus_nhet-1]++;
					clus_nhet=0;
				}				
				temp_concentration_HetR+=(*it).Get_HetR();
			}
			
			if((it-filament.begin())==0)
			{
				temp_nbhet=0;
				std::vector<CCell>::iterator tit=it;
				while ( !(*tit).Get_vegetative() && tit<filament.end() )
				{
					temp_nbhet++;
					tit++;
					if((tit-filament.end())+1==0) fullH=1;
				}
				if(temp_nbhet<BORDER_SIZE-1) temp_Bhist[temp_nbhet]++;
				else temp_Bhist[BORDER_SIZE-1]++;
			}
			else if(((it-filament.end())+1==0)&&(fullH==0))
			{
				if (clus_nhet!=0)
				{
					temp_chet[clus_nhet-1]++;
					clus_nhet=0;
				}	
				temp_nbhet=0;
				std::vector<CCell>::iterator tit=it;
				while ( !(*tit).Get_vegetative() && tit>=filament.begin() )
				{
					temp_nbhet++;
					tit--;
				}
				if(temp_nbhet<BORDER_SIZE-1) temp_Bhist[temp_nbhet]++;
				else temp_Bhist[BORDER_SIZE-1]++;
			}
		}

		temp_ncells=filament.size();
		n_cells_med[indice_nt]+=temp_ncells;
		n_cells_disp[indice_nt]+=gsl_pow_2(temp_ncells);			
		n_clushet_med[indice_nt]+=temp_nchet;
		n_clushet_disp[indice_nt]+=gsl_pow_2(temp_nchet);						
		
		if(temp_nhet<temp_ncells) concentration_HetR_filament_med[indice_nt]+=temp_concentration_HetR/(temp_ncells-temp_nhet); //Concentración media de HetR total en las células vegetativas
		n_het_med[indice_nt]+=temp_nhet;
		n_het_disp[indice_nt]+=gsl_pow_2(temp_nhet);
		n_veg_med[indice_nt]+=(temp_ncells-temp_nhet);
		n_veg_disp[indice_nt]+=gsl_pow_2(temp_ncells-temp_nhet);
		n_het_boundaries_med[indice_nt]+=temp_nbhet;
		n_het_boundaries_disp[indice_nt]+=gsl_pow_2(temp_nbhet);
		
		n_cells_batch_med[kSim*((int)(Ttot/AT)+1)+indice_nt]+=temp_ncells;			
		n_het_batch_med[kSim*((int)(Ttot/AT+1))+indice_nt]+=temp_nhet;
		n_het_boundaries_batch_med[kSim*((int)(Ttot/AT)+1)+indice_nt]+=temp_nbhet;		
		n_clushet_batch_med[kSim*((int)(Ttot/AT)+1)+indice_nt]+=temp_nchet;		
			
		p_clushet_med[indice_nt]+=(double)temp_nchet/temp_ncells;
		p_clushet_disp[indice_nt]+=gsl_pow_2((double)temp_nchet/temp_ncells);
		p_clushet_batch_med[kSim*((int)(Ttot/AT)+1)+indice_nt]+=(double)temp_nchet/temp_ncells;		
		
		p_het_med[indice_nt]+=(double)temp_nhet/temp_ncells;
		p_het_disp[indice_nt]+=gsl_pow_2((double)temp_nhet/temp_ncells);
		p_het_batch_med[kSim*((int)(Ttot/AT+1))+indice_nt]+=(double)temp_nhet/temp_ncells;	

		for(int i=0;i<BORDER_SIZE;i++)
		{
			n_het_bound_hist_med[i+indice_nt*BORDER_SIZE]+=temp_Bhist[i];      //Suma en cada realización el valor para obtener la media
			n_het_bound_hist_disp[i+indice_nt*BORDER_SIZE]+=gsl_pow_2(temp_Bhist[i]); //Suma en cada realización el cuadrado para obtener la dispersión			
			n_het_bound_hist_batch_med[NBch*(i+indice_nt*BORDER_SIZE)+kSim]+=temp_Bhist[i];      //Suma en cada realización el valor para obtener la dispersión por batches
			
			clushet_size_med[i+indice_nt*BORDER_SIZE]+=temp_chet[i];
			clushet_size_disp[i+indice_nt*BORDER_SIZE]+=gsl_pow_2(temp_chet[i]);			
			clushet_size_batch_med[NBch*(i+indice_nt*BORDER_SIZE)+kSim]+=temp_chet[i];
		}

		for(int i=0;i<HISTOGRAM_SIZE;i++)
		{
			histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]+=temp_histogram[i];      //Suma en cada realización el valor para obtener la media
			histogram_veg_interval_disp[i+indice_nt*HISTOGRAM_SIZE]+=gsl_pow_2(temp_histogram[i]); //Suma en cada realización el cuadrado para obtener la dispersión						
			histogram_veg_interval_batch_med[NBch*(i+indice_nt*HISTOGRAM_SIZE)+kSim]+=temp_histogram[i];      //Suma en cada realización el valor para obtener la media
		}
	}
	
	void Batch_dispersion()
	{
		for(int indice_nt=0; indice_nt<=Ttot/AT; indice_nt++)
		{
			for(int k=0;k<NBch;k++)
			{	
				n_cells_batch_disp[indice_nt]+=gsl_pow_2(n_cells_batch_med[k*((int)(Ttot/AT+1))+indice_nt]);
				n_het_batch_disp[indice_nt]+=gsl_pow_2(n_het_batch_med[k*((int)(Ttot/AT)+1)+indice_nt]);
				n_het_boundaries_batch_disp[indice_nt]+=gsl_pow_2(n_het_boundaries_batch_med[k*((int)(Ttot/AT)+1)+indice_nt]);				
				n_clushet_batch_disp[indice_nt]+=gsl_pow_2(n_clushet_batch_med[k*((int)(Ttot/AT+1))+indice_nt]);
				//~ std::cout<<indice_nt<<", "<<k<<"\nMed: "<<n_clushet_batch_med[k*((int)(Ttot/AT+1))+indice_nt]<<"\t"<<"Disp: "<<n_clushet_batch_disp[indice_nt]<<"\n";	
				
				p_het_batch_disp[indice_nt]+=gsl_pow_2(p_het_batch_med[k*((int)(Ttot/AT)+1)+indice_nt]);				
				p_clushet_batch_disp[indice_nt]+=gsl_pow_2(p_clushet_batch_med[k*((int)(Ttot/AT+1))+indice_nt]);		
				
				for(int i=0;i<HISTOGRAM_SIZE;i++)
				{
					histogram_veg_interval_batch_disp[i+indice_nt*HISTOGRAM_SIZE]+=gsl_pow_2(histogram_veg_interval_batch_med[NBch*(i+indice_nt*HISTOGRAM_SIZE)+k]);
					if(i<BORDER_SIZE)
					{
						n_het_bound_hist_batch_disp[i+indice_nt*BORDER_SIZE]+=gsl_pow_2(n_het_bound_hist_batch_med[NBch*(i+indice_nt*BORDER_SIZE)+k]);
						clushet_size_batch_disp[i+indice_nt*BORDER_SIZE]+=gsl_pow_2(clushet_size_batch_med[NBch*(i+indice_nt*BORDER_SIZE)+k]);
					}
				}					
			}
		}
	}
	
	void Print_histogram(std::string Mpath)
	{
		//~ for(int indice_nt=0; indice_nt<=Ttot/AT; indice_nt++)
		//~ for(int indice_nt=12/AT; indice_nt<=Ttot/AT; indice_nt+=12/AT)
		for(int indice_nt=30/AT; indice_nt<=60/AT; indice_nt+=20/AT)
		{
			int total=0, totalB=0, totalCH=0;
			std::string Hpath=Mpath+"/Histograms/histogram_med_t="+DoubleToString(indice_nt*AT,2)+".dat";
			std::string HBpath=Mpath+"/Histograms/Bhistogram_med_t="+DoubleToString(indice_nt*AT,2)+".dat";
			std::string CHpath=Mpath+"/Histograms/CHhistogram_med_t="+DoubleToString(indice_nt*AT,2)+".dat";
			std::ofstream f1(Hpath);
			std::ofstream f2(HBpath);
			std::ofstream f3(CHpath);

			for(int i=0;i<HISTOGRAM_SIZE;i++)
			{
				total+=histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]; //obtienen el numero total de intervalos para normalizar
				if(i<BORDER_SIZE)
				{
					totalB+=n_het_bound_hist_med[i+indice_nt*BORDER_SIZE];
					totalCH+=clushet_size_med[i+indice_nt*BORDER_SIZE];					
				}			
			}

			for(int i=0;i<HISTOGRAM_SIZE;i++)
			{
				if(total>0)	f1 <<i<<"\t\t"<<histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]/total<<"\t\t"<<sqrt((histogram_veg_interval_disp[i+indice_nt*HISTOGRAM_SIZE]-gsl_pow_2(histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]/MED)*MED)/(MED-1))*MED/total<<"\t\t"<<sqrt((histogram_veg_interval_batch_disp[i+indice_nt*HISTOGRAM_SIZE]-gsl_pow_2(histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]/NBch)*NBch)/(NBch-1))*NBch/total<<std::endl;
				else f1 <<i<<"\t\t"<<0 <<"\t\t"<<0 <<"\t\t"<<0<<std::endl;

				if(i<BORDER_SIZE)
				{
					if(totalB>0)	f2<<i<<"\t\t"<<n_het_bound_hist_med[i+indice_nt*BORDER_SIZE]/totalB<<"\t\t"<<sqrt((n_het_bound_hist_disp[i+indice_nt*BORDER_SIZE]-gsl_pow_2(n_het_bound_hist_med[i+indice_nt*BORDER_SIZE]/MED)*MED)/(MED-1))*MED/total<<"\t\t"<<sqrt((n_het_bound_hist_batch_disp[i+indice_nt*BORDER_SIZE]-gsl_pow_2(n_het_bound_hist_med[i+indice_nt*BORDER_SIZE]/NBch)*NBch)/(NBch-1))*NBch/total<<std::endl;
					else f2 <<i<<"\t\t"<<0 <<"\t\t"<<0 <<"\t\t"<<0<<std::endl;
					if(totalCH>0)	f3<<i+1<<"\t\t"<<clushet_size_med[i+indice_nt*BORDER_SIZE]/totalCH<<"\t\t"<<sqrt((clushet_size_disp[i+indice_nt*BORDER_SIZE]-gsl_pow_2(clushet_size_med[i+indice_nt*BORDER_SIZE]/MED)*MED)/(MED-1))*MED/total<<"\t\t"<<sqrt((clushet_size_batch_disp[i+indice_nt*BORDER_SIZE]-gsl_pow_2(clushet_size_med[i+indice_nt*BORDER_SIZE]/NBch)*NBch)/(NBch-1))*NBch/total<<std::endl;
					else f3 <<i<<"\t\t"<<0 <<"\t\t"<<0 <<"\t\t"<<0<<std::endl;
				}
			}

			f1<<std::endl<<"#-------------------------------------"<<std::endl<<"#Intervals \t frecuency \t deviation \t Batch deviation"<<std::endl;
			f1<<"#Averaged="<<MED<<" times."<<std::endl;
			f1<<"#Batches="<<NBch<<" times."<<std::endl;
			f1.close();

			f2<<std::endl<<"#-------------------------------------"<<std::endl<<"#Intervals \t frecuency \t deviation \t Batch deviation"<<std::endl;
			f2<<"#Averaged="<<MED<<" times."<<std::endl;
			f2<<"#Batches="<<NBch<<" times."<<std::endl;
			f2.close();
			
			f3<<std::endl<<"#-------------------------------------"<<std::endl<<"#Cluster Size \t frecuency \t deviation \t Batch deviation"<<std::endl;
			f3<<"#Averaged="<<MED<<" times."<<std::endl;
			f3<<"#Batches="<<NBch<<" times."<<std::endl;
			f3.close();
		}
	}

	void Print_data_histogram(std::string Mpath)
	{
		std::string mean=Mpath+"/histogram_mean.dat";
		std::string variance=Mpath+"/histogram_variance.dat";
		std::string skewness=Mpath+"/histogram_skewness.dat";
		std::string median=Mpath+"/histogram_median.dat";
		std::string evens=Mpath+"/even_intervals.dat";
		std::string contiguous=Mpath+"/histogram_contiguous.dat";

		std::ofstream f1(mean);
		std::ofstream f2(variance);
		std::ofstream f3(skewness);
		std::ofstream f4(median);
		std::ofstream f5(evens);
		std::ofstream f6(contiguous);
		
		for(int indice_nt=0; indice_nt<=Ttot/AT; indice_nt++)
		{
			double mean=0, median=0, variance=0, skewness=0, deviation_mean=0, Batch_deviation_mean=0;
			int total=0, total_evens=0, total_intervals_less_cero=0, temp=0;

			for(int i=0;i<HISTOGRAM_SIZE-1;i++) //Va recorriendo el histograma para cada instante de tiempo menos la última posición donde se guarda el acumulativo (por eso -1)
			{					//Para los datos estadísticos se deja fuera el acumulado del histograma
				total+=histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE];
				mean+=i*histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]; //El error de la media es la raiz cuadrada de la suma de los errores al cuadrado (hay que dividir luego por el total)
				deviation_mean+=i*i*(histogram_veg_interval_disp[i+indice_nt*HISTOGRAM_SIZE]-gsl_pow_2(histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]/MED)*MED)/(MED-1);
				Batch_deviation_mean+=i*i*(histogram_veg_interval_batch_disp[i+indice_nt*HISTOGRAM_SIZE]-gsl_pow_2(histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]/NBch)*NBch)/(NBch-1);
				variance+=i*i*histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE];
				skewness+=i*i*i*histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE];
				if(i!=0)
				{
					total_intervals_less_cero+=histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE];
					if(i%2==0) total_evens+=histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE]; //Guarda el número de intervalos pares
				}
			}

			for(int i=0;temp<total/2. && i<HISTOGRAM_SIZE-1;i++)   //recorre el histograma hasta que temp es menor que el total
			{
				temp+=histogram_veg_interval_med[i+indice_nt*HISTOGRAM_SIZE];
				median++;
			}

			if(total>0)
			{
				double An_cont=sqrt((histogram_veg_interval_disp[indice_nt*HISTOGRAM_SIZE]-gsl_pow_2(histogram_veg_interval_med[indice_nt*HISTOGRAM_SIZE]/MED)*MED)/(MED-1))*MED/total;
				double Batch_An_cont=sqrt((histogram_veg_interval_batch_disp[indice_nt*HISTOGRAM_SIZE]-gsl_pow_2(histogram_veg_interval_med[indice_nt*HISTOGRAM_SIZE]/NBch)*NBch)/(NBch-1))*NBch/total;
				
				f1<<indice_nt*AT<<"\t\t"<<mean/total<<"\t\t"<<sqrt(deviation_mean)*MED/total<<"\t\t"<<sqrt(deviation_mean)*NBch/total<<std::endl;
				f2<<indice_nt*AT<<"\t\t"<<variance/total-gsl_pow_2(mean/total)<<std::endl;
				f3<<indice_nt*AT<<"\t\t"<<(skewness/total-3*mean/total*variance/total+2*gsl_pow_3(mean/total))/gsl_pow_3(sqrt(variance/total-gsl_pow_2(mean/total)))<<std::endl;
				if(total%2==0) f4<<indice_nt*AT<<"\t\t"<<median-0.5<<std::endl;
				else f4<<indice_nt*AT<<"\t\t"<<median-1<<std::endl;
				f5<<indice_nt*AT<<"\t\t"<<total_evens/total_intervals_less_cero<<std::endl;
				
				f6<<indice_nt*AT<<"\t\t"<<histogram_veg_interval_med[indice_nt*HISTOGRAM_SIZE]/total <<"\t\t"<<An_cont<<"\t\t"<<Batch_An_cont<<std::endl;
			}
			else
			{
				f1<<indice_nt*AT<<"\t\t"<<0<<"\t\t"<<0<<"\t\t"<<0<<std::endl;
				f2<<indice_nt*AT<<"\t\t"<<0<<std::endl;
				f3<<indice_nt*AT<<"\t\t"<<0<<std::endl;
				f4<<indice_nt*AT<<"\t\t"<<0<<std::endl;
				f5<<indice_nt*AT<<"\t\t"<<0<<std::endl;
				f6<<indice_nt*AT<<"\t\t"<<0<<"\t\t"<<0<<"\t\t"<<0<<std::endl;
			}
		}

		f1<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time \t mean \t deviation \t Batch deviation"<<std::endl;
		f1<<"#Averaged="<<MED<<" times."<<std::endl;
		f1<<"#Batches="<<NBch<<" times."<<std::endl;		
		f1.close();

		f2<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time \t variance"<<std::endl;
		f2<<"#Averaged="<<MED<<" times."<<std::endl;
		f2.close();

		f3<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time \t skewness"<<std::endl;
		f3<<"#Averaged="<<MED<<" times."<<std::endl;
		f3.close();

		f4<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time \t median"<<std::endl;
		f4<<"#Averaged="<<MED<<" times."<<std::endl;
		f4.close();

		f5<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time \t even intervals"<<std::endl;
		f5<<"#Averaged="<<MED<<" times."<<std::endl;
		f5.close();

		f6<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time \t contiguous heterocyst (%) \t deviation \t Batch deviation"<<std::endl;
		f6<<"#Averaged="<<MED<<" times."<<std::endl;
		f6<<"#Batches="<<NBch<<" times."<<std::endl;		
		f6.close();
	}

	void Print_cells(std::string Mpath)
	{
		std::string percent=Mpath+"/per_cent_het.dat";
		std::string ncells=Mpath+"/n_cells_total.dat";
		std::string nhet=Mpath+"/n_heterocists.dat";
		std::string nhetbord=Mpath+"/n_heterocists_boundaries.dat";
		std::string nveg=Mpath+"/n_vegetatives.dat";
		std::string concveg=Mpath+"/concentration_HetR.dat";
		std::string nclusters=Mpath+"/histogram_nclusters.dat";
		std::string pclusters=Mpath+"/per_cent_clusters.dat";
		std::ofstream f1(percent);
		std::ofstream f2(ncells);
		std::ofstream f3(nhet);
		std::ofstream f4(nveg);
		std::ofstream f5(nhetbord);
		std::ofstream f6(concveg);
		std::ofstream f7(nclusters);
		std::ofstream f8(pclusters);

		for(int i=0;i<=Ttot/AT;i++)
		{
			double An_het=sqrt((n_het_disp[i]-gsl_pow_2(n_het_med[i]/MED)*MED)/(MED-1));
			double Batch_An_het=sqrt((n_het_batch_disp[i]-gsl_pow_2(n_het_med[i]/NBch)*NBch)/(NBch-1))*NBch/MED;
			double An_cells=sqrt((n_cells_disp[i]-gsl_pow_2(n_cells_med[i]/MED)*MED)/(MED-1));
			double Batch_An_cells=sqrt((n_cells_batch_disp[i]-gsl_pow_2(n_cells_med[i]/NBch)*NBch)/(NBch-1))*NBch/MED;
			double An_clus=sqrt((n_clushet_disp[i]-gsl_pow_2(n_clushet_med[i]/MED)*MED)/(MED-1));
			double Batch_An_clus=sqrt((n_clushet_batch_disp[i]-gsl_pow_2(n_clushet_med[i]/NBch)*NBch)/(NBch-1))*NBch/MED;
			//~ double An_percH=n_het_med[i]/n_cells_med[i]*sqrt(gsl_pow_2(An_het/n_het_med[i])+gsl_pow_2(An_cells/n_cells_med[i]))*100
			//~ double Batch_An_percH=n_het_med[i]/n_cells_med[i]*sqrt(gsl_pow_2(Batch_An_het/n_het_med[i])+gsl_pow_2(Batch_An_cells/n_cells_med[i]))*100
			double An_percH=sqrt((p_het_disp[i]-gsl_pow_2(p_het_med[i]/MED)*MED)/(MED-1))*100;
			double Batch_An_percH=sqrt((p_het_batch_disp[i]-gsl_pow_2(p_het_med[i]/NBch)*NBch)/(NBch-1))*NBch/MED*100;
			double An_perC=sqrt((p_clushet_disp[i]-gsl_pow_2(p_clushet_med[i]/MED)*MED)/(MED-1))*100;
			double Batch_An_perC=sqrt((p_clushet_batch_disp[i]-gsl_pow_2(p_clushet_med[i]/NBch)*NBch)/(NBch-1))*NBch/MED*100;

			if(n_het_med[i]!=0)
			{
				//~ f1<<i*AT<<"\t\t"<<n_het_med[i]/n_cells_med[i]*100<<"\t\t"<<An_percH<<"\t\t"<<Batch_An_percH<<std::endl;
				f1<<i*AT<<"\t\t"<<p_het_med[i]/MED*100<<"\t\t"<<An_percH<<"\t\t"<<Batch_An_percH<<std::endl;
				f7<<i*AT<<"\t\t"<<n_clushet_med[i]/MED<<"\t\t"<<An_clus<<"\t\t"<<Batch_An_clus<<std::endl;
				f8<<i*AT<<"\t\t"<<p_clushet_med[i]/MED*100<<"\t\t"<<An_perC<<"\t\t"<<Batch_An_perC<<std::endl;
			}
			else
			{
				f1<<i*AT<<"\t\t"<<0<<"\t\t"<<0<<"\t\t"<<0<<std::endl;
				if (n_het_med[i]==n_cells_med[i])
				{
					f7<<i*AT<<"\t\t"<<1<<"\t\t"<<0<<"\t\t"<<0<<std::endl;
					f8<<i*AT<<"\t\t"<<1<<"\t\t"<<0<<"\t\t"<<0<<std::endl;
				}
				else
				{
					f7<<i*AT<<"\t\t"<<0<<"\t\t"<<0<<"\t\t"<<0<<std::endl;
					f8<<i*AT<<"\t\t"<<0<<"\t\t"<<0<<"\t\t"<<0<<std::endl;
				}
			}

			f2<<i*AT<<"\t\t"<<n_cells_med[i]/MED<<"\t\t"<<An_cells<<"\t\t"<<Batch_An_cells<<std::endl;
			f3<<i*AT<<"\t\t"<<n_het_med[i]/MED<<"\t\t"<<An_het<<"\t\t"<<Batch_An_het<<std::endl;
			f4<<i*AT<<"\t\t"<<n_veg_med[i]/MED<<"\t\t"<<sqrt((n_veg_disp[i]-gsl_pow_2(n_veg_med[i]/MED)*MED)/(MED-1))<<"\t\t"<<sqrt(gsl_pow_2(Batch_An_het)+gsl_pow_2(Batch_An_cells))<<std::endl;
			f5<<i*AT<<"\t\t"<<n_het_boundaries_med[i]/MED<<"\t\t"<<sqrt((n_het_boundaries_disp[i]-gsl_pow_2(n_het_boundaries_med[i]/MED)*MED)/(MED-1))<<"\t\t"<<sqrt((n_het_boundaries_batch_disp[i]-gsl_pow_2(n_het_boundaries_med[i]/NBch)*NBch)/(NBch-1))*NBch/MED<<std::endl;
			f6<<i*AT<<"\t\t"<<concentration_HetR_filament_med[i]/MED<<std::endl;
		}

		f1<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t n_heterocysts(%) \t deviation \t Batch deviation"<<std::endl;
		f1<<"#Averaged="<<MED<<" times."<<std::endl;
		f1<<"#Batches="<<NBch<<" times."<<std::endl;		
		f1.close();

		f2<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t Averaged n. of cells \t deviation \t Batch deviation"<<std::endl;
		f2<<"#Averaged="<<MED<<" times."<<std::endl;
		f2<<"#Batches="<<NBch<<" times."<<std::endl;		
		f2.close();

		f3<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t Averaged n. of heterocists \t deviation \t Batch deviation"<<std::endl;
		f3<<"#Averaged="<<MED<<" times."<<std::endl;
		f3<<"#Batches="<<NBch<<" times."<<std::endl;		
		f3.close();

		f4<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t Averaged n. of vegetative cells \t deviation \t Batch deviation"<<std::endl;
		f4<<"#Averaged="<<MED<<" times."<<std::endl;
		f4<<"#Batches="<<NBch<<" times."<<std::endl;		
		f4.close();

		f5<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t Averaged n. of heterocists at the boundaries \t deviation \t Batch deviation"<<std::endl;
		f5<<"#Averaged="<<MED<<" times."<<std::endl;
		f5<<"#Batches="<<NBch<<" times."<<std::endl;		
		f5.close();

		f6<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t Averaged total concentration of HetR in vegetative cells \t Batch deviation"<<std::endl;
		f6<<"#Averaged="<<MED<<" times."<<std::endl;
		f6<<"#Batches="<<NBch<<" times."<<std::endl;		
		f6.close();

		f7<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t Averaged number of heterocysts clusters \t deviation \t Batch deviation"<<std::endl;
		f7<<"#Averaged="<<MED<<" times."<<std::endl;
		f7<<"#Batches="<<NBch<<" times."<<std::endl;		
		f7.close();
		
		f8<<std::endl<<"#-------------------------------------"<<std::endl<<"#Time(h.) \t n_clusters(%) \t deviation \t Batch deviation"<<std::endl;
		f8<<"#Averaged="<<MED<<" times."<<std::endl;
		f8<<"#Batches="<<NBch<<" times."<<std::endl;		
		f8.close();
	}

	//~ void GFP_position(std::vector<CCell> &filament, int indice_nt)
	//~ {
		//~ int BegTag=0, position=0;
		//~ int PosFil[filament.size()];
		//~ int TPosFil[filament.size()];
		//~ int i;
		//~ for(i=0;i<filament.size();i++)
		//~ {
			//~ PosFil[i]=PMAX+1;
			//~ TPosFil[i]=TPMAX+1;
		//~ }

		//~ for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)   //El "for" lo avanza antes. Está bien
		//~ {
			//~ position++;
			//~ if(BegTag==1) PosFil[it-filament.begin()]=position;

			//~ if(!(*it).Get_vegetative())
			//~ {
				//~ BegTag=1;
				//~ position=0;
				//~ PosFil[it-filament.begin()]=0;
				//~ if((it-filament.begin())!=0)
				//~ {
					//~ position++;
					//~ while((it-filament.begin()-position)!=0 && position<=PosFil[it-filament.begin()-position])
					//~ {
						//~ PosFil[it-filament.begin()-position]=position;
						//~ position++;
					//~ }
					//~ position=0;
				//~ }
			//~ }
		//~ }

		//~ if(PosFil[0]==0) //Terminal Tag
		//~ {
			//~ TPosFil[0]=0;
			//~ if(TPMAX>-1) PosFil[0]=PMAX+1; //Only Internal
			//~ i=1;
			//~ while(PosFil[i]!=0 && i<TPMAX+1)
			//~ {
				//~ if (TPMAX>0 && i<=PosFil[i]) PosFil[i]=PMAX+1; //Only Internal
				//~ TPosFil[i]=i;
				//~ i++;
			//~ }
		//~ }

		//~ if(PosFil[filament.size()-1]==0) //Terminal Tag
		//~ {
			//~ TPosFil[filament.size()-1]=0;
			//~ if(TPMAX>-1) PosFil[filament.size()-1]=PMAX+1; //Only Internal
			//~ i=1;
			//~ while(PosFil[filament.size()-1-i]!=0 && i<TPMAX+1)
			//~ {
				//~ if (TPMAX>0 && i<=PosFil[i]) PosFil[filament.size()-1-i]=PMAX+1; //Only Internal
				//~ TPosFil[filament.size()-1-i]=i;
				//~ i++;
			//~ }
		//~ }

		//~ // char control[256];
		//~ // strcpy(control, name_directory);
		//~ // strcat(control, "/controlPTAG.dat");
		//~ // std::ofstream f9(control);
		//~ //
		//~ // for(int j=0; j<filament.size(); j++){f9<<PosFil[j]<<"\t";}
		//~ // f9<<std::endl;
		//~ // for(int j=0; j<filament.size(); j++){f9<<TPosFil[j]<<"\t";}
		//~ // f9<<std::endl;
		//~ //
		//~ // f9<<std::endl<<"#-------------------------------------"<<std::endl;
		//~ // for(int j=0; j<filament.size(); j++){f9<<"I"<<j<<"\t";}
		//~ // f9<<std::endl;
		//~ // for(int j=0; j<filament.size(); j++){f9<<"T"<<j<<"\t";}
		//~ // f9<<std::endl;
		//~ // f9.close();

		//~ for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)   //El "for" lo avanza antes. Está bien
		//~ {
			//~ if (PosFil[it-filament.begin()]<PMAX+1)
			//~ {
				//~ GFP_ConPos_n[indice_nt][PosFil[it-filament.begin()]]++;
				//~ GFP_ConPos_med[indice_nt][PosFil[it-filament.begin()]]+=(*it).Get_GFP();
				//~ GFP_ConPos_disp[indice_nt][PosFil[it-filament.begin()]]+=gsl_pow_2((*it).Get_GFP());
				//~ HetR_ConPos_n[indice_nt][PosFil[it-filament.begin()]]++;
				//~ HetR_ConPos_med[indice_nt][PosFil[it-filament.begin()]]+=(*it).Get_HetR();
				//~ HetR_ConPos_disp[indice_nt][PosFil[it-filament.begin()]]+=gsl_pow_2((*it).Get_HetR());
			//~ }

			//~ if (TPosFil[it-filament.begin()]<TPMAX+1)
			//~ {
				//~ TGFP_ConPos_n[indice_nt][TPosFil[it-filament.begin()]]++;
				//~ TGFP_ConPos_med[indice_nt][TPosFil[it-filament.begin()]]+=(*it).Get_GFP();
				//~ TGFP_ConPos_disp[indice_nt][TPosFil[it-filament.begin()]]+=gsl_pow_2((*it).Get_GFP());
				//~ THetR_ConPos_n[indice_nt][TPosFil[it-filament.begin()]]++;
				//~ THetR_ConPos_med[indice_nt][TPosFil[it-filament.begin()]]+=(*it).Get_HetR();
				//~ THetR_ConPos_disp[indice_nt][TPosFil[it-filament.begin()]]+=gsl_pow_2((*it).Get_HetR());
			//~ }
		//~ }
	//~ }

	//~ void Print_GFP_position(std::string Mpath, char name_directory[128])
	//~ {
		//~ std::string time_means_GFP=Mpath+"/time_means_GFP.dat";
		//~ std::string time_means_HetR=Mpath+"/time_means_HetR.dat";
		//~ std::string CompSimGFP=Mpath+"/CompSimGFP.dat";
		//~ std::string CompSimHetR=Mpath+"/CompSimHetR.dat";
		//~ std::ofstream f1(time_means_GFP);
		//~ std::ofstream f2(CompSimGFP);
		//~ std::ofstream f3(time_means_HetR);
		//~ std::ofstream f4(CompSimHetR);

		//~ double allmean_GFP[PMAX+1], alldisp_GFP[PMAX+1], alln_GFP[PMAX+1], allmean_HetR[PMAX+1], alldisp_HetR[PMAX+1], alln_HetR[PMAX+1];
		//~ double normGFP=0, normHetR=0, Ntotal=0;
		//~ for(int j=0; j<PMAX+1; j++)
		//~ {
			//~ allmean_GFP[j]=0;
			//~ alldisp_GFP[j]=0;
			//~ alln_GFP[j]=0;
			//~ allmean_HetR[j]=0;
			//~ alldisp_HetR[j]=0;
			//~ alln_HetR[j]=0;
		//~ }

		//~ char terminal_means_GFP[256];
		//~ char terminal_means_HetR[256];
		//~ char TCompSimGFP[256];
		//~ char TCompSimHetR[256];
		//~ if(TPMAX>-1)
		//~ {
			//~ strcpy(terminal_means_GFP, name_directory);
			//~ strcat(terminal_means_GFP, "/terminal_means_GFP.dat");
			//~ strcpy(terminal_means_HetR, name_directory);
			//~ strcat(terminal_means_HetR, "/terminal_means_HetR.dat");
			//~ strcpy(TCompSimGFP, name_directory);
			//~ strcat(TCompSimGFP, "/TCompSimGFP.dat");
			//~ strcpy(TCompSimHetR, name_directory);
			//~ strcat(TCompSimHetR, "/TCompSimHetR.dat");
		//~ }
		//~ std::ofstream f5(terminal_means_GFP);
		//~ std::ofstream f6(TCompSimGFP);
		//~ std::ofstream f7(terminal_means_HetR);
		//~ std::ofstream f8(TCompSimHetR);

		//~ double Tallmean_GFP[TPMAX+1], Talldisp_GFP[TPMAX+1], Talln_GFP[TPMAX+1], Tallmean_HetR[TPMAX+1], Talldisp_HetR[TPMAX+1], Talln_HetR[TPMAX+1];
		//~ double normTGFP=0, normTHetR=0, TNtotal=0;
		//~ for(int j=0; j<TPMAX+1; j++)
		//~ {
			//~ Tallmean_GFP[j]=0;
			//~ Talldisp_GFP[j]=0;
			//~ Talln_GFP[j]=0;
			//~ Tallmean_HetR[j]=0;
			//~ Talldisp_HetR[j]=0;
			//~ Talln_HetR[j]=0;
		//~ }

		//~ for(int i=0; i<=Ttot/AT; i++)
		//~ {
			//~ for(int j=0; j<PMAX+1; j++)
			//~ {
				//~ allmean_GFP[j]+=GFP_ConPos_med[i][j];
				//~ alldisp_GFP[j]+=GFP_ConPos_disp[i][j];
				//~ alln_GFP[j]+=GFP_ConPos_n[i][j];
				//~ normGFP+=GFP_ConPos_med[i][j];
				//~ allmean_HetR[j]+=HetR_ConPos_med[i][j];
				//~ alldisp_HetR[j]+=HetR_ConPos_disp[i][j];
				//~ alln_HetR[j]+=HetR_ConPos_n[i][j];
				//~ normHetR+=HetR_ConPos_med[i][j];
				//~ Ntotal+=GFP_ConPos_n[i][j];
			//~ }
			//~ if(TPMAX>-1)
			//~ {
				//~ for(int j=0; j<PMAX+1; j++)
				//~ {
					//~ Tallmean_GFP[j]+=TGFP_ConPos_med[i][j];
					//~ Talldisp_GFP[j]+=TGFP_ConPos_disp[i][j];
					//~ Talln_GFP[j]+=TGFP_ConPos_med[i][j];
					//~ normTGFP+=Tallmean_GFP[j]*Talln_GFP[j];
					//~ Tallmean_HetR[j]+=THetR_ConPos_med[i][j];
					//~ Talldisp_HetR[j]+=THetR_ConPos_disp[i][j];
					//~ Talln_HetR[j]+=THetR_ConPos_n[i][j];
					//~ normTHetR+=THetR_ConPos_med[i][j];
					//~ TNtotal+=TGFP_ConPos_n[i][j];
				//~ }
			//~ }
		//~ }
		//~ normGFP/=Ntotal;
		//~ normHetR/=Ntotal;
		//~ normTGFP/=TNtotal;
		//~ normTHetR/=TNtotal;

		//~ for(int i=0; i<=Ttot/AT; i++)
		//~ {
			//~ f1<<i*AT;
			//~ for(int j=0; j<PMAX+1; j++){f1<<"\t"<<(GFP_ConPos_med[i][j]/GFP_ConPos_n[i][j])/normGFP<<"\t"<<sqrt((GFP_ConPos_disp[i][j]-gsl_pow_2(GFP_ConPos_med[i][j]/GFP_ConPos_n[i][j])*GFP_ConPos_n[i][j])/(GFP_ConPos_n[i][j]-1))/normGFP<<"\t"<<GFP_ConPos_n[i][j];}
			//~ f1<<std::endl;

			//~ f3<<i*AT;
			//~ for(int j=0; j<PMAX+1; j++){f3<<"\t"<<(HetR_ConPos_med[i][j]/HetR_ConPos_n[i][j])/normHetR<<"\t"<<sqrt((HetR_ConPos_disp[i][j]-gsl_pow_2(HetR_ConPos_med[i][j]/HetR_ConPos_n[i][j])*HetR_ConPos_n[i][j])/(HetR_ConPos_n[i][j]-1))/normHetR<<"\t"<<HetR_ConPos_n[i][j];}
			//~ f3<<std::endl;

			//~ if(TPMAX>-1)
			//~ {
				//~ f5<<i*AT;
				//~ for(int j=0; j<TPMAX+1; j++){f5<<"\t"<<(TGFP_ConPos_med[i][j]/TGFP_ConPos_n[i][j])/normTGFP<<"\t"<<sqrt((TGFP_ConPos_disp[i][j]-gsl_pow_2(TGFP_ConPos_med[i][j]/TGFP_ConPos_n[i][j])*TGFP_ConPos_n[i][j])/(TGFP_ConPos_n[i][j]-1))/normTGFP<<"\t"<<TGFP_ConPos_n[i][j];}
				//~ f5<<std::endl;

				//~ f7<<i*AT;
				//~ for(int j=0; j<TPMAX+1; j++){f7<<"\t"<<(THetR_ConPos_med[i][j]/THetR_ConPos_n[i][j])/normTHetR<<"\t"<<sqrt((THetR_ConPos_disp[i][j]-gsl_pow_2(THetR_ConPos_med[i][j]/THetR_ConPos_n[i][j])*THetR_ConPos_n[i][j])/(THetR_ConPos_n[i][j]-1))/normTHetR<<"\t"<<THetR_ConPos_n[i][j];}
				//~ f7<<std::endl;
			//~ }
		//~ }

		//~ f1<<std::endl<<"#-------------------------------------"<<std::endl;
		//~ f1<<"#Time";
		//~ for(int j=0; j<PMAX+1; j++){f1<<"\t"<<"Mean P"<<j<<"\t"<<"Deviation P"<<j<<"\t"<<"#Cases P"<<j;}
		//~ f1<<std::endl;
		//~ f1<<"#Simulations="<<MED<<std::endl;
		//~ f1.close();

		//~ f3<<std::endl<<"#-------------------------------------"<<std::endl;
		//~ f3<<"#Time";
		//~ for(int j=0; j<PMAX+1; j++){f3<<"\t"<<"Mean P"<<j<<"\t"<<"Deviation P"<<j<<"\t"<<"#Cases P"<<j;}
		//~ f3<<std::endl;
		//~ f3<<"#Simulations="<<MED<<std::endl;
		//~ f3.close();

		//~ if(TPMAX>-1)
		//~ {
			//~ f5<<std::endl<<"#-------------------------------------"<<std::endl;
			//~ f5<<"#Time";
			//~ for(int j=0; j<TPMAX+1; j++){f5<<"\t"<<"Mean P"<<j<<"\t"<<"Deviation P"<<j<<"\t"<<"#Cases P"<<j;}
			//~ f5<<std::endl;
			//~ f5<<"#Simulations="<<MED<<std::endl;
			//~ f5.close();

			//~ f7<<std::endl<<"#-------------------------------------"<<std::endl;
			//~ f7<<"#Time";
			//~ for(int j=0; j<TPMAX+1; j++){f7<<"\t"<<"Mean P"<<j<<"\t"<<"Deviation P"<<j<<"\t"<<"#Cases P"<<j;}
			//~ f7<<std::endl;
			//~ f7<<"#Simulations="<<MED<<std::endl;
			//~ f7.close();
		//~ }

		//~ for(int j=0; j<PMAX+1; j++)
		//~ {
			//~ f2<<j<<"\t\t"<<(allmean_GFP[j]/alln_GFP[j])/normGFP<<"\t\t"<<sqrt((alldisp_GFP[j]-gsl_pow_2(allmean_GFP[j]/alln_GFP[j])*alln_GFP[j])/(alln_GFP[j]-1))/normGFP<<"\t\t"<<sqrt(((alldisp_GFP[j]-gsl_pow_2(allmean_GFP[j]/alln_GFP[j])*alln_GFP[j])/(alln_GFP[j]-1))/alln_GFP[j])/normGFP<<std::endl;

			//~ f4<<j<<"\t\t"<<(allmean_HetR[j]/alln_HetR[j])/normHetR<<"\t\t"<<sqrt((alldisp_HetR[j]-gsl_pow_2(allmean_HetR[j]/alln_HetR[j])*alln_HetR[j])/(alln_HetR[j]-1))/normHetR<<"\t\t"<<sqrt(((alldisp_HetR[j]-gsl_pow_2(allmean_HetR[j]/alln_HetR[j])*alln_HetR[j])/(alln_HetR[j]-1))/alln_HetR[j])/normHetR<<std::endl;
		//~ }

		//~ f2<<std::endl<<"#-------------------------------------"<<std::endl<<"#Position \t Mean \t STDeviation \t STError"<<std::endl;
		//~ for(int j=0; j<PMAX+1; j++){f2<<"#P"<<j<<"\t"<<alln_GFP[j]<<std::endl;}
		//~ f2<<"#Simulations="<<MED<<std::endl;
		//~ f2.close();

		//~ f4<<std::endl<<"#-------------------------------------"<<std::endl<<"#Position \t Mean \t STDeviation \t STError"<<std::endl;
		//~ for(int j=0; j<PMAX+1; j++){f4<<"#P"<<j<<"\t"<<alln_HetR[j]<<std::endl;}
		//~ f4<<"#Simulations="<<MED<<std::endl;
		//~ f4.close();

		//~ if(TPMAX>-1)
		//~ {
			//~ for(int j=0; j<TPMAX+1; j++)
			//~ {
				//~ f6<<j<<"\t\t"<<(Tallmean_GFP[j]/Talln_GFP[j])/normTGFP<<"\t\t"<<sqrt((Talldisp_GFP[j]-gsl_pow_2(Tallmean_GFP[j]/Talln_GFP[j])*Talln_GFP[j])/(Talln_GFP[j]-1))/normTGFP<<"\t\t"<<sqrt(((Talldisp_GFP[j]-gsl_pow_2(Tallmean_GFP[j]/Talln_GFP[j])*Talln_GFP[j])/(Talln_GFP[j]-1))/Talln_GFP[j])/normTGFP<<std::endl;

				//~ f8<<j<<"\t\t"<<(Tallmean_HetR[j]/Talln_HetR[j])/normTHetR<<"\t\t"<<sqrt((Talldisp_HetR[j]-gsl_pow_2(Tallmean_HetR[j]/Talln_HetR[j])*Talln_HetR[j])/(Talln_HetR[j]-1))/normTHetR<<"\t\t"<<sqrt(((Talldisp_HetR[j]-gsl_pow_2(Tallmean_HetR[j]/Talln_HetR[j])*Talln_HetR[j])/(Talln_HetR[j]-1))/Talln_HetR[j])/normTHetR<<std::endl;
			//~ }

			//~ f6<<std::endl<<"#-------------------------------------"<<std::endl<<"#Position \t Mean \t STDeviation \t STError"<<std::endl;
			//~ for(int j=0; j<TPMAX+1; j++){f6<<"#P"<<j<<"\t"<<Talln_GFP[j]<<std::endl;}
			//~ f6<<"#Simulations="<<MED<<std::endl;
			//~ f6.close();

			//~ f8<<std::endl<<"#-------------------------------------"<<std::endl<<"#Position \t Mean \t STDeviation \t STError"<<std::endl;
			//~ for(int j=0; j<TPMAX+1; j++){f8<<"#P"<<j<<"\t"<<Talln_HetR[j]<<std::endl;}
			//~ f8<<"#Simulations="<<MED<<std::endl;
			//~ f8.close();
		//~ }
	//~ }
};

void Initialize_filament(std::vector<CCell> &filament, int N, const std::vector<double>& Param)
{
	filament.clear();
	for (int i=0; i<N; i++)
	{
		CCell tempcell(Param);
		filament.push_back(tempcell);
	}
}

std::vector<CCell>::iterator Insert_cell(std::vector<CCell> &filament, std::vector<CCell>::iterator it)
{
   return ++filament.insert(it, *it); 			//Duplica la celula y devuelve un iterator a la segunda copia insertada
}

void Advance_At(std::vector<CCell> &filament)
{
	double size2;
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		if((*it).Get_protoheterocist())
		{
			//std::cout<<"Cell"<<it-filament.begin()<<"\n";
			(*it).Advance_At_protoheterocist(filament, it-filament.begin() );	//Avanca un paso de tiempo todas las variables (tamaño y concentraciones)

			if((*it).Get_PROTO_T()>T_MEAN_PROTO)  //Si supera el nivel acumulativo intrínseco max_level durante T_MEAN convierte a protoheterociste
			{
				(*it).Transform_heterocist();
				//std::cout<<it-filament.begin()<<", ";
			}
		}
		else if((*it).Get_vegetative())
		{
			//std::cout<<"Cell"<<it-filament.begin()<<"\n";
			(*it).Advance_At_vegetative(filament, it-filament.begin() );	//Avanca un paso de tiempo todas las variables (tamaño y concentraciones)

			if((*it).Get_size() > (*it).GeT_MEAN_size())
			{
				double OHetRA=(*it).Get_cumulative_HetR();
				size2=(*it).Set_size_cell1();
				(*it).Set_cell(OHetRA);        //Establece threshholds de la primera célula hija
				it=Insert_cell(filament, it); //Duplica la célula si ha alcanzado el tamaño máximo y coloca el iterator en la siguiente
				//std::cout<<"Cell "<<it-filament.begin()<<", "<<(*it).Get_HetR()<<", "<<"Cell "<<it-1-filament.begin()<<", "<<(*(it-1)).Get_HetR()<<"\n";
				(*it).Set_cell(OHetRA);        //Establece threshholds de la segunda célula hija
				(*it).Set_size_cell2(size2);
			}

			if(((*it).Get_DifD_T()>T_MIN)&&((*it).Get_cumulative_HetR()>(*it).GeT_MEAN_level()*T_MEAN))  //Si supera el nivel acumulativo intrínseco max_level durante T_MEAN convierte a heterociste
			{
				if (T_MEAN_PROTO>At)
					(*it).Transform_protoheterocist();
				else
					(*it).Transform_heterocist();
			}
		}
		else
		{
			(*it).Advance_At_heterocist(filament, it-filament.begin());
			//if(it-filament.begin()==56)  std::cout<< (*it).Get_HetR()<<", ";
		}
	}
}

void Cout_filament(std::vector<CCell> &filament) //WTF ¿?
{
	//std::cout<<std::endl;
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		std::cout <<(*it).Get_HetR()<<", ";
	}
	std::cout<<std::endl;
}

void Print_filament_HetR(std::vector<CCell> &filament, double nt, std::string Mpath)
{
	std::string Fpath=Mpath+"/Profiles_HetR/Profile_HetR_t="+DoubleToString(nt*At,2)+".dat";
	std::ofstream fs(Fpath);

	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		fs <<it-filament.begin()<<"\t\t"<<(*it).Get_HetR()<<std::endl;
	}
	fs<<std::endl<<"#-------------------------------------"<<std::endl<<"# i \t HetR(i)"<<std::endl;
	fs.close();
}

//~ void Print_filament_GFP(std::vector<CCell> &filament, double nt, std::string Mpath)
//~ {
	//~ std::string Fpath=Mpath+"/Profiles_GFP/Profile_GFP_t="+DoubleToString(nt*At,2)+".dat";
	//~ std::ofstream fs(Fpath);

	//~ for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	//~ {
		//~ fs <<it-filament.begin()<<"\t\t"<<(*it).Get_GFP()<<std::endl;
	//~ }
	//~ fs<<std::endl<<"#-------------------------------------"<<std::endl<<"# i \t GFP(i)"<<std::endl;
	//~ fs.close();
//~ }

void Print_filament_PatA(std::vector<CCell> &filament, double nt, std::string Mpath)
{
	std::string Fpath=Mpath+"/Profiles_PatA/Profile_PatA_t="+DoubleToString(nt*At,2)+".dat";
	std::ofstream fs(Fpath);

	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		fs <<it-filament.begin()<<"\t\t"<<(*it).Get_PatA()<<std::endl;
	}
	fs<<std::endl<<"#-------------------------------------"<<std::endl<<"# i \t PatA(i)"<<std::endl;
	fs.close();
}

void Print_filament_PatS(std::vector<CCell> &filament, double nt, std::string Mpath)
{
	std::string Fpath=Mpath+"/Profiles_PatS/Profile_PatS_t="+DoubleToString(nt*At,2)+".dat";
	std::ofstream fs(Fpath);

	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		fs <<it-filament.begin()<<"\t\t"<<(*it).Get_PatS()<<std::endl;
	}
	fs<<std::endl<<"#-------------------------------------"<<std::endl<<"# i \t PatS(i)"<<std::endl;
	fs.close();
}

void Print_filament_HetN(std::vector<CCell> &filament, double nt, std::string Mpath)
{
	std::string Fpath=Mpath+"/Profiles_HetN/Profile_HetN_t="+DoubleToString(nt*At,2)+".dat";
	std::ofstream fs(Fpath);

	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		fs <<it-filament.begin()<<"\t\t"<<(*it).Get_HetN()<<std::endl;
	}
	fs<<std::endl<<"#-------------------------------------"<<std::endl<<"# i \t HetN(i)"<<std::endl;
	fs.close();
}

void Print_filament_Inhb(std::vector<CCell> &filament, double nt, std::string Mpath)
{
	std::string Fpath=Mpath+"/Profiles_Inhb/Profile_Inhb_t="+DoubleToString(nt*At,2)+".dat";
	std::ofstream fs(Fpath);

	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		fs <<it-filament.begin()<<"\t\t"<<(*it).Get_Inhb()<<std::endl;
	}
	fs<<std::endl<<"#-------------------------------------"<<std::endl<<"# i \t Inhb(i)"<<std::endl;
	fs.close();
}

void Print_filament_fixN(std::vector<CCell> &filament, double nt, std::string Mpath)
{
	std::string Fpath=Mpath+"/Profiles_fixN/Profile_fixN_t="+DoubleToString(nt*At,2)+".dat";
	std::ofstream fs(Fpath);

	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		fs <<it-filament.begin()<<"\t\t"<<(*it).Get_fixN()<<std::endl;
	}
	fs<<std::endl<<"#-------------------------------------"<<std::endl<<"# i \t fixN(i)"<<std::endl;
	fs.close();
}

void DrawText(SDL_Surface* screen, TTF_Font* font, const char* text, Sint16 posx, Sint16 posy) //Función para dibujar el texto
{
	SDL_Color color={255,255,255};
	SDL_Rect font_position={posx,posy};
	SDL_Surface *text_surface;

	text_surface=TTF_RenderText_Solid(font, text, color);

	SDL_BlitSurface(text_surface, NULL, screen, &font_position);
	SDL_FreeSurface(text_surface);
}

void Draw_filament(std::vector<CCell> &filament, double nt, SDL_Surface* screen, std::string Mpath)
{
	int cell_size=12;
	double MAX_LEVEL_HetR=100;
	double MIN_LEVEL_HetR=15;
	double MAX_LEVEL_Inhb=20;
	double MIN_LEVEL_Inhb=5;
	double MAX_LEVEL_fixN=150;
	double MIN_LEVEL_fixN=20;

	std::string Fpath=Mpath+"/Figuras/Filament_t="+DoubleToString(nt*At,2)+".bmp";
	std::string Title="Time "+DoubleToString(nt*At,2)+" h.";

	char name_filament[512], font_text[64];
	strcpy(name_filament, Fpath.c_str());
	strcpy(font_text, Title.c_str());

	SDL_WM_SetCaption(name_filament, 0);
	SDL_Event event;
	SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 0, 0, 0)); //Rellena la pantalla. (0, 0, 0) de negro, (255, 255, 255) de blanco

	//Escribe tiempo pantalla.
	TTF_Init();
	TTF_Font* font_time=TTF_OpenFont("Arial.ttf", 40);
	DrawText(screen, font_time, font_text, WINDOW_WIDTH/2-100, 5); ////Llama a la función para dibujar el texto y lo coloca en posx posy
	TTF_CloseFont(font_time);
	TTF_Quit();

	//Escribe HetR en la pantalla.
	TTF_Init();
	TTF_Font* font_protein=TTF_OpenFont("Arial.ttf", 32);
	DrawText(screen, font_protein, "HetR", 10, 40); ////Llama a la función para dibujar el texto y lo coloca en posx posy

	double radio, level, position=20, position_row=90;

	//Dibuja las vegetativas en la pantalla con la concentración de HetR
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		if(position>=WINDOW_WIDTH-40) //Coloca la posición en la siguiente linea si alcanza el final
		{	position_row+=50;   position=20;    }

		if((*it).Get_vegetative())
		{
			level=255*((*it).Get_HetR()-MIN_LEVEL_HetR)/MAX_LEVEL_HetR; 
			if(level>255) level=255;
			else if (level<0) level=0;
			
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			//filledEllipseRGBA(screen, position, position_row, radio+1, cell_size+1, 255, 255, 255, 255);  //Color de la membrana de la célula
			//filledEllipseRGBA(screen, position, position_row, radio, cell_size+1, 29, 91, 56, 255);   	//Color de la membrana de la célula como en el esquema
			if ((*it).Get_protoheterocist())
				filledEllipseRGBA(screen, position, position_row, radio+1.5, cell_size+1.5, 255, 0, 0, 150);   		//Color de la membrana de la célula rojo
			else
				filledEllipseRGBA(screen, position, position_row, radio, cell_size, 255, 255, 255, 150);   		//Color de la membrana de la célula blanco
			filledEllipseRGBA(screen, position, position_row, radio-1, cell_size-1, 0, 0, 0, 255);   		//Color de fondo de la célula (negro)
			filledEllipseRGBA(screen, position, position_row, radio-1, cell_size-1, 0, 192 , 0, level);  	//Color de la concentración de proteína
			position+=radio;
		}
		else
		{
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			position+=radio;
		}
	}

	position=20; position_row=90; //Reposiciona al inicio

	//Dibuja los heterocistes en la pantalla con la concentración de HetR
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		if(position>=WINDOW_WIDTH-40) //Coloca la posición en la siguiente linea si alcanza el final
		{   position_row+=50;   position=20;    }

		if((*it).Get_vegetative())
		{
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			position+=radio;
		}
		else
		{
			level=255; if(level>255) level=255;
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			filledEllipseRGBA(screen, position, position_row, radio+0.25+1, 14.5+1, 255 ,255, 255, 255);
			filledEllipseRGBA(screen, position, position_row, radio+0.25-1, 14.5-1, 0 ,0, 0, 255);	   //Color de fondo de la célula (negro)
			filledEllipseRGBA(screen, position, position_row, radio+0.25-1, 14.5-1, 0, 192 , 0, level);  //Color de la concentración de proteína
			position+=radio;
		}
	}
	double position_row_init=position_row;
	
	//Escribe Inhb en la pantalla.
	position=20;    //Establece de nuevo la posicion del nuevo filamento
	position_row=position_row+50;
	DrawText(screen, font_protein, "Inhb", 10, position_row); ////Llama a la función para dibujar el texto y lo coloca en posx posy
	position_row=position_row+50;
	position_row_init=position_row;

	//Dibuja las vegetativas en la pantalla con la concentración de Inhb
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		if(position>=WINDOW_WIDTH-40) //Coloca la posición en la siguiente linea si alcanza el final
		{   position_row+=50;	position=20;  	}

		if((*it).Get_vegetative())
		{
			level=255*((*it).Get_Inhb()-MIN_LEVEL_Inhb)/MAX_LEVEL_Inhb; 
			if(level>255) level=255;
			else if (level<0) level=0;
				
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			//filledEllipseRGBA(screen, position, position_row, radio+1, cell_size+1, 255, 255, 255, 255);   //Color de la membrana de la célula
			//filledEllipseRGBA(screen, position, position_row, radio, cell_size+1, 29, 91, 56, 255);   //Color de la membrana de la célula como en el esquema
			if ((*it).Get_protoheterocist())
				filledEllipseRGBA(screen, position, position_row, radio+1.5, cell_size+1.5, 255, 0, 0, 150);   		//Color de la membrana de la célula rojo
			else
				filledEllipseRGBA(screen, position, position_row, radio, cell_size, 255, 255, 255, 150);   		//Color de la membrana de la célula blanco
			filledEllipseRGBA(screen, position, position_row, radio-1, cell_size-1, 0, 0, 0, 255);   //Color de fondo de la célula (negro)
			filledEllipseRGBA(screen, position, position_row, radio-1, cell_size-1, 150, 1 , 200, level);  //Color de la concentración de proteína
			position+=radio;
		}
		else
		{
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			position+=radio;
		}
	}

	position=20; position_row=position_row_init; //Reposiciona al inicio

	//Dibuja los heterocistes en la pantalla con la concentración de Inhb
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		if(position>=WINDOW_WIDTH-40) //Coloca la posición en la siguiente linea si alcanza el final
		{   position_row+=50;   position=20;    }

		if((*it).Get_vegetative())
		{
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			position+=radio;
		}
		else
		{
			level=255*((*it).Get_Inhb()-MIN_LEVEL_Inhb)/MAX_LEVEL_Inhb; 
			if(level>255) level=255;
			else if (level<0) level=0;
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			filledEllipseRGBA(screen, position, position_row, radio+0.25+1, 14.5+1, 255 ,255, 255, 255);
			filledEllipseRGBA(screen, position, position_row, radio+0.25-1, 14.5-1, 0 ,0, 0, 255);	   //Color de fondo de la célula (negro)
			filledEllipseRGBA(screen, position, position_row, radio+0.25-1, 14.5-1, 150, 1 , 200, level);  //Color de la concentración de proteína
			position+=radio;
		}
	}

	//Escribe fixN en la pantalla.
	position=20;    //Establece de nuevo la posicion del nuevo filamento
	position_row=position_row+50;
	DrawText(screen, font_protein, "fixed N", 10, position_row); ////Llama a la función para dibujar el texto y lo coloca en posx posy
	position_row=position_row+50;
	position_row_init=position_row;

	//Dibuja las vegetativas en la pantalla con la concentración de fixN
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		if(position>=WINDOW_WIDTH-40) //Coloca la posición en la siguiente linea si alcanza el final
		{   position_row+=50;   position=20;    }

		if((*it).Get_vegetative())
		{
			level=255*((*it).Get_fixN()-MIN_LEVEL_fixN)/MAX_LEVEL_fixN; 
			if(level>255) level=255;
			else if (level<0) level=0;
			
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			//filledEllipseRGBA(screen, position, position_row, radio+1, cell_size+1, 255, 255, 255, 255);   //Color de la membrana de la célula
			//filledEllipseRGBA(screen, position, position_row, radio, cell_size+1, 29, 91, 56, 255);   //Color de la membrana de la célula como en el esquema
			if ((*it).Get_protoheterocist())
				filledEllipseRGBA(screen, position, position_row, radio+1.5, cell_size+1.5, 255, 0, 0, 150);   		//Color de la membrana de la célula rojo
			else
				filledEllipseRGBA(screen, position, position_row, radio, cell_size, 255, 255, 255, 150);   		//Color de la membrana de la célula blanco
			filledEllipseRGBA(screen, position, position_row, radio-1, cell_size-1, 0, 0, 0, 255);   //Color de fondo de la célula (negro)
			filledEllipseRGBA(screen, position, position_row, radio-1, cell_size-1, 0, 255 , 255, level);  //Color de la concentración de proteína
		position+=radio;
		}
		else
		{
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			position+=radio;
		}
	}

	position=20; position_row=position_row_init; //Reposiciona al inicio

	//Dibuja los heterocistes en la pantalla con la concentración de fixN
	for (std::vector<CCell>::iterator it=filament.begin(); it<filament.end(); it++)
	{
		if(position>=WINDOW_WIDTH-40) //Coloca la posición en la siguiente linea si alcanza el final
		{   position_row+=50;   position=20;    }

		if((*it).Get_vegetative())
		{
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			position+=radio;
		}
		else
		{
			level=255*((*it).Get_fixN()-MIN_LEVEL_fixN)/MAX_LEVEL_fixN; 
			if(level>255) level=255;
			else if (level<0) level=0;
			radio=14/MAX_SIZE*(*it).Get_size();
			position+=radio;
			filledEllipseRGBA(screen, position, position_row, radio+0.25+1, 14.5+1, 255 ,255, 255, 255);
			filledEllipseRGBA(screen, position, position_row, radio+0.25-1, 14.5-1, 0 ,0, 0, 255);	   //Color de fondo de la célula (negro)
			filledEllipseRGBA(screen, position, position_row, radio+0.25-1, 14.5-1, 0, 255 , 255, level);  //Color de la concentración de proteína
			position+=radio;
		}
	}

	TTF_CloseFont(font_protein);
	TTF_Quit();

	SDL_SaveBMP(screen, name_filament);
	SDL_Flip(screen);
	// SDL_Delay( 5000 );
}

void build_directory(std::string path)
{
	std::string cmd = "mkdir "+ path;
	system (cmd.c_str());
}

std::string create_directories(std::string SParam,std::string Mutant)
{
	std::string Mpath="Data_"+SParam+"_"+Mutant;
	build_directory(Mpath);

	std::string path=Mpath+"/Figuras";
	build_directory(path);

	path=Mpath+"/Profiles_HetR";
	build_directory(path);

	//~ path=Mpath+"/Profiles_GFP";
	//~ build_directory(path);

	path=Mpath+"/Profiles_PatA";
	build_directory(path);

	path=Mpath+"/Profiles_PatS";
	build_directory(path);

	path=Mpath+"/Profiles_HetN";
	build_directory(path);

	path=Mpath+"/Profiles_Inhb";
	build_directory(path);

	path=Mpath+"/Profiles_fixN";
	build_directory(path);

	path=Mpath+"/Histograms";
	build_directory(path);

	return(Mpath);
}

std::vector<double> read_IParam(std::string SParam)
{
	std::vector<double> Param;
	if(SParam.at(0) == 'N')
	{
		// if first character is an N then delete the N
		// but only delete the N if it's the first character
		SParam.erase(SParam.begin()+0);
	}
	std::string filename="Parameters/BestParam"+SParam+".txt";
	std::ifstream in(filename);
	double r2;
	while (in >> r2)//cuando llega al final del archivo, esta condición devuelve False, y se sale del bucle
	{
		Param.push_back(r2);
	}
	in.close();//cierra archivo
	return Param;
}

int main(int argc, char *argv[])
{
	std::stringstream sSParam;
	std::string SParam;

	sSParam << argv[1];
	sSParam >> SParam;

	r=gsl_rng_alloc (gsl_rng_mt19937); 		//Inicializamos el generador de números aleatorios
	gsl_rng_set(r, SEMILLA);				 //Inicializamos con nuestra semilla

	int n_At=(int)(AT/At);  				//Numero de pasos temporales para representar datos

	std::vector<std::string> Mutants;
	Mutants.push_back("WT");
	Mutants.push_back("DPS");
	Mutants.push_back("DHN");
	Mutants.push_back("DPA");
	Mutants.push_back("DPAHN");
	Mutants.push_back("DPX");
	Mutants.push_back("DPAPS");
	Mutants.push_back("DHF");
	//~ Mutants.push_back("FDPS");
	//~ //Noise variations
	//~ Mutants.push_back("WTLCN");
	//~ Mutants.push_back("WTLRN");
	//~ Mutants.push_back("WTLN");
	//~ Mutants.push_back("WTHCN");
	//~ Mutants.push_back("WTHRN");
	//~ Mutants.push_back("WTHN");
	//~ //patS/patX variations
	//~ Mutants.push_back("DPS2");
	//~ Mutants.push_back("DPX2");	
	//~ //Difusion variations
	//~ Mutants.push_back("DPAFD");
	//~ Mutants.push_back("DPAND");	
	//~ Mutants.push_back("DPAHNFD");
	//~ Mutants.push_back("DPAHNND");

	char name_directory[128];
	std::string Mpath;

	SDL_Init( SDL_INIT_VIDEO );
	SDL_Surface* screen=SDL_SetVideoMode( WINDOW_WIDTH, WINDOW_HEIGHT, 32, SDL_HWSURFACE | SDL_ANYFORMAT ); //*/
	std::cout.setf(std::ios::fixed);   // Hace que no se almacenen las instrucciones del std::cout en el bufer sino que las imprime
	std::cout.setf(std::ios::unitbuf); //instantaneamente por pantalla cuando le llegan

	double ControlBatch[NBch];
	for (int k=0;k<NBch;k++)
		ControlBatch[k]=0;
		
	for(auto Mutant: Mutants)
	{
		std::vector<double> Param=read_IParam(SParam);

		//~ SDL_Init( SDL_INIT_VIDEO );
		//~ SDL_Surface* screen=SDL_SetVideoMode( WINDOW_WIDTH, WINDOW_HEIGHT, 32, SDL_HWSURFACE | SDL_ANYFORMAT ); //*/
		//~ std::cout.setf(std::ios::fixed);   // Hace que no se almacenen las instrucciones del std::cout en el bufer sino que las imprime
		//~ std::cout.setf(std::ios::unitbuf); //instantaneamente por pantalla cuando le llegan
		
		//Basic Mutants
		if(Mutant=="WT")
		{
			Mpath=create_directories(SParam,"WT");
		}
		else if(Mutant=="DPS")
		{
			Mpath=create_directories(SParam,"patSdel");
			Param[10]=PatX*Param[10];
		}
		else if(Mutant=="DHN")
		{
			Mpath=create_directories(SParam,"hetNdel");
			Param[13]=0;
		}
		else if(Mutant=="DPA")
		{
			Mpath=create_directories(SParam,"patAdel");
			Param[6]=0;
		}
		else if(Mutant=="DPAHN")
		{
			Mpath=create_directories(SParam,"patAhetNdel");
			Param[6]=0;
			Param[13]=0;
		}
		else if(Mutant=="DPX")
		{
			Mpath=create_directories(SParam,"patXdel");
			Param[10]=(1-PatX)*Param[10];
		}
		else if(Mutant=="DPAPS")
		{
			Mpath=create_directories(SParam,"patApatSdel");
			Param[6]=0;
			Param[10]=PatX*Param[10];
		}
		else if(Mutant=="DHF")
		{
			Mpath=create_directories(SParam,"hetFdel");
			Param[9]=0;
		}
		//Noise study
		else if(Mutant=="WTLCN")
		{
			Mpath=create_directories(SParam,"WTLCN");
			Param[24]=0.1*Param[24];
		}
		else if(Mutant=="WTLRN")
		{
			Mpath=create_directories(SParam,"WTLRN");
			Param[25]=0.1*Param[25];
		}
		else if(Mutant=="WTLN")
		{
			Mpath=create_directories(SParam,"WTLN");
			//~ Param[24]=0.1*Param[24];
			Param[25]=0.1*Param[25];
		}
		else if(Mutant=="WTHCN")
		{
			Mpath=create_directories(SParam,"WTHCN");
			Param[24]=10*Param[24];
		}
		else if(Mutant=="WTHRN")
		{
			Mpath=create_directories(SParam,"WTHRN");
			Param[25]=4*Param[25];
		}
		else if(Mutant=="WTHN")
		{
			Mpath=create_directories(SParam,"WTHN");
			Param[24]=10*Param[24];
			Param[25]=4*Param[25];
		}
		//patS variants
		else if(Mutant=="FDPS")
		{
			Mpath=create_directories(SParam,"patXpatSdel");
			Param[10]=0;
		}
		else if(Mutant=="DPS2")
		{
			Mpath=create_directories(SParam,"patSdel2");
			Param[10]=2*PatX*Param[10];
		}
		else if(Mutant=="DPX2")
		{
			Mpath=create_directories(SParam,"patXdel2");
			Param[10]=(1-PatX*2)*Param[10];
		}
		//Difusion variants
		else if(Mutant=="DPAFD")
		{
			Mpath=create_directories(SParam,"patAdel FullDifussion");
			Param[6]=0;
			Param[23]=1;
		}
		else if(Mutant=="DPAND")
		{
			Mpath=create_directories(SParam,"patAdel NoDifussion");
			Param[6]=0;
			Param[23]=0;
		}
		else if(Mutant=="DPAHNFD")
		{
			Mpath=create_directories(SParam,"patAhetNdel FullDifussion");
			Param[6]=0;
			Param[13]=0;
			Param[23]=1;
		}
		else if(Mutant=="DPAHNND")
		{
			Mpath=create_directories(SParam,"patAhetNdel NoDifussion");
			Param[6]=0;
			Param[13]=0;
			Param[23]=0;
		}
		strcpy(name_directory,Mpath.c_str());

		CHistogram histogram;
		int NumThr;
		if(MED<MaxThr) NumThr=MED;
		else NumThr=MaxThr;
		#pragma omp parallel for num_threads(NumThr)
		for(int k=1;k<=MED;++k)
		{
			//~ std::cout<<"\t"<<k<<"\t"<<k%NBch<<"\n"<<std::endl;
			std::vector<CCell> filament;
			T_MEAN_PROTO=Param[26];
			Initialize_filament(filament, INITIAL_N_CELLS, Param);
			for(long int nt=0;nt<=Ttot/At;nt++)
			{
				if(nt%n_At==0)
				{
					if(k==1)
					{
						Draw_filament(filament,nt,screen,Mpath);
						Print_filament_HetR(filament,nt,Mpath);
						//~ Print_filament_GFP(filament,nt,Mpath);
						Print_filament_PatA(filament,nt,Mpath);
						Print_filament_PatS(filament,nt,Mpath);
						Print_filament_HetN(filament,nt,Mpath);
						Print_filament_Inhb(filament,nt,Mpath);
						Print_filament_fixN(filament,nt,Mpath);
					}
					
					histogram.Sum_histogram(filament,nt/n_At,k%NBch);
					//~ histogram.GFP_position(filament,nt/n_At);
				}
				Advance_At(filament);
			}
			//std::cout<<"variable="<<rho_s<<", Med="<<k<<", nº total células="<<filament.size()<<"; ";
			filament.clear();
			if(Mutant=="WT") ControlBatch[k%NBch]++;
		}
		histogram.Batch_dispersion();
		//~ histogram.Print_GFP_position(Mpath,name_directory);
		histogram.Print_data_histogram(Mpath);
		histogram.Print_cells(Mpath);
		histogram.Print_histogram(Mpath);
		histogram.CHistogram::~CHistogram();  //Delete object histogram
		//~ SDL_Quit();
	}
	SDL_Quit();	
	for (int k=0;k<NBch;k++)
		std::cout<<"Batch["<<k<<"] size: "<<ControlBatch[k]<<std::endl;
	return 0;
}
