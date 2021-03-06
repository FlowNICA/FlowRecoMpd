#ifndef UTILITY_H

#define UTILITY_H

#define _MAX_TRACKS 5000
#define _N_ARM 2 // 3 IS FULL DETECTOR
#define _N_HARM 2
#define _N_SORTS 4
#define _N_MODULES_TOTAL 90
#define _N_METHOD 2 // 0 - TPC , 1 - ZDC
#define _N_QCOMP 2

#include <TMath.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TApplication.h>
#include "SpecFuncMathMore.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>

R__LOAD_LIBRARY(libMathMore.so);

using TMath::Abs;
using TMath::Cos;
using TMath::Sin;
using TMath::ATan2;
using TMath::Sqrt;
using TMath::LocMin;
using TMath::Pi;

using std::cout;
using std::endl;

const Float_t Cut_Pt_Min = 0.;
const Float_t Cut_Eta_Min = 0.7;
const Float_t Cut_Eta_Max = 1.5;
const Float_t Cut_Eta_Gap = 0.05;
const Int_t Cut_No_Of_hits_min = 32;

const int Ndim = 3;

const float centralityBinsFlow[] = {0.,10.,40.,80};
const int NcentralityBinsFlow = 3;

//const float centralityBinsRes[] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100};
//const int NcentralityBinsRes = 10;
const float centralityBinsRes[] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100};
const int NcentralityBinsRes = 20;

const double ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const int NptBins = 12;

const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const int NetaBins = 14;

const float rapidityBins[] = {-2., -1.8, -1.6, -1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.};
const int NrapidityBins = 20;

const TString arm_names[_N_ARM] = {TString("R"),TString("L")};
const TString sorts_of_particles[4] = {TString("all sorts") , TString("pions (211)") , TString("protons (2212)") , TString("kaons (321)")};
const TString methods_names[_N_METHOD] = {TString("TPC") , TString("FHCal")};

class FlowParticle
{
public:

    double Eta;
    double Pt;
    double Phi;
    double Rapidity;

    FlowParticle();

    FlowParticle(double Eta, double Pt, double Phi, double Rapidity);
};

class EPParticle
{
public:

    double Eta;
    double Pt;
    double Phi;

    EPParticle();

    EPParticle(double Eta, double Pt, double Phi);
};

Double_t ResEventPlane(Double_t chi, Int_t harm); //harm = 1 or 2 for our case

Double_t Chi(Double_t res, Int_t harm); //harm = 1 or 2 for our case

Double_t* GetAngles();

Float_t Unfold(Float_t phiEP_mc, Float_t psi_N_FULL, Int_t harm);

#endif
