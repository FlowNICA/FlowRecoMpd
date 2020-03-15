
#include <cstdio>

#include <TROOT.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TF1.h>
#include <TFile.h>

#include "utility.C"

/*const int _N_HARM = 2;
const int _N_METHOD = 2;
const int _N_SORTS = 4;

const float centralityBinsFlow[] = {10.,20.,40.,50};
const int NcentralityBinsFlow = 3;

const float centralityBinsRes[] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100};
const int NcentralityBinsRes = 20;

const float ptBins[] = {0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const int NptBins = 12;

const float etaBins[] = {-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const int NetaBins = 14;

const float rapidityBins[] = {-2., -1.8, -1.6, -1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.};
const int NrapidityBins = 20;

const TString sorts_of_particles[4] = {TString("all sorts") , TString("pions (#pi^{+})") , TString("protons (p)") , TString("kaons (K^{+})")};

const TString methods_names[_N_METHOD] = {TString("TPC") , TString("FHCal")};*/

using TMath::Sqrt;

//R__ADD_LIBRARY_PATH("/lustre/nyx/hades/user/parfenov/mpd_new/real-flow") // if needed
//R__LOAD_LIBRARY(MpdCalculator_C.so)
//R__LOAD_LIBRARY(utility_C.so)

TGraph *Convert2Graph(TH1D *const& prof)
{
  std::vector<Double_t> vX, vY;
  double y=0.;
  vX = {5.,15.,25.,35.,45.,55.,65.,75.,85.,95.};
  for (int i=0; i<prof->GetNbinsX()/2; i++)
  {
    y = prof->GetBinContent(2*i+1);
    if (y <= 0.) continue;
    vY.push_back(y);
  }
  TGraph *graph = new TGraph(vY.size(),&vX[0],&vY[0]);
  return graph;
}

void get_res(TString inFileName , TString outFileName)
{
	
	//gSystem->Load("./utility_C.so");
	//gSystem->Load("./MpdCalculator_C.so");
		
	TFile *inFile = new TFile(inFileName.Data());
	TFile *outFile = new TFile(outFileName.Data(),"RECREATE");
	
	TProfile *p_Res2Psi_vs_b[_N_HARM][_N_HARM][_N_METHOD];
	TH1D *p_Res[_N_HARM][_N_HARM][_N_METHOD];
	TF1 *resolution_fit[_N_HARM][_N_HARM][_N_METHOD];
	
	//TH1F *h_pt_PID_efficiency_before[NetaBins];
	//TH1F *h_pt_PID_efficiency_after_sorts[NetaBins][_N_SORTS];
	//TH1F *h_pt_PID_eff[NetaBins][_N_SORTS];
	//TF1  *f_pt_PID_eff[NetaBins][_N_SORTS];
	
	TH1F *h_pt_reco_before[NcentralityBinsRes][_N_SORTS][NetaBins];
	TH1F *h_pt_true_before[NcentralityBinsRes][_N_SORTS][NetaBins];
	TH1F *h_pt_reco_after[NcentralityBinsRes][_N_SORTS][NetaBins];
	TH1F *h_pt_true_after[NcentralityBinsRes][_N_SORTS][NetaBins];
	
	TH1F *h_pt_eff_before[NcentralityBinsRes][_N_SORTS][NetaBins];
	TH1F *h_pt_eff_after[NcentralityBinsRes][_N_SORTS][NetaBins];
	TF1  *f_pt_eff_before[NcentralityBinsRes][_N_SORTS][NetaBins];
	TF1  *f_pt_eff_after[NcentralityBinsRes][_N_SORTS][NetaBins];
	
	TH1F *h_pt_reco_final[NcentralityBinsFlow][_N_SORTS][NetaBins];
	TH1F *h_pt_true_final[NcentralityBinsFlow][_N_SORTS][NetaBins];
	TH1F *h_pt_eff_after_final[NcentralityBinsFlow][_N_SORTS][NetaBins];
	TF1  *f_pt_eff_after_final[NcentralityBinsFlow][_N_SORTS][NetaBins];
	
	TH1F *h_eta_reco_before[NcentralityBinsRes][_N_SORTS][NptBins];
	TH1F *h_eta_true_before[NcentralityBinsRes][_N_SORTS][NptBins];
	TH1F *h_eta_reco_after[NcentralityBinsRes][_N_SORTS][NptBins];
	TH1F *h_eta_true_after[NcentralityBinsRes][_N_SORTS][NptBins];
	
	TH1F *h_eta_eff_before[NcentralityBinsRes][_N_SORTS][NptBins];
	TH1F *h_eta_eff_after[NcentralityBinsRes][_N_SORTS][NptBins];
	TF1  *f_eta_eff_before[NcentralityBinsRes][_N_SORTS][NptBins];
	TF1  *f_eta_eff_after[NcentralityBinsRes][_N_SORTS][NptBins];
	
	TH1F *h_eta_reco_final[NcentralityBinsFlow][_N_SORTS][NptBins];
	TH1F *h_eta_true_final[NcentralityBinsFlow][_N_SORTS][NptBins];
	TH1F *h_eta_eff_after_final[NcentralityBinsFlow][_N_SORTS][NptBins];
	TF1  *f_eta_eff_after_final[NcentralityBinsFlow][_N_SORTS][NptBins];
	
	char name[400];
	char title[400];
	
	//for (Int_t eta_bin=0; eta_bin<NetaBins; eta_bin++){
		//h_pt_PID_efficiency_before[eta_bin] = (TH1F*) inFile->Get(Form("h_pt_PID_efficiency_before%i",eta_bin));
		
		//for (int sort = 0; sort < _N_SORTS; ++sort){
			//sprintf(title,"p_{T} efficiency of PID cuts for %s, %.2f < #eta < %.2f ", sorts_of_particles[sort].Data(), etaBins[eta_bin], etaBins[eta_bin+1]);
			//h_pt_PID_efficiency_after_sorts[eta_bin][sort] = (TH1F*) inFile->Get(Form("h_pt_PID_efficiency_after_sorts%i%i",eta_bin,sort));
			
			//h_pt_PID_eff[eta_bin][sort] = new TH1F(Form("h_pt_PID_eff%i%i",eta_bin,sort),title,NptBins,ptBins); h_pt_PID_eff[eta_bin][sort] -> Sumw2();
			//f_pt_PID_eff[eta_bin][sort] = new TF1(Form("f_pt_PID_eff%i%i",eta_bin,sort),"pol8",0.2,2.);
		//}
	//}
	
	/*for (int centralityBin = 0; centralityBin < NcentralityBinsFlow; ++centralityBin)
	{
		for (int sort = 0; sort < _N_SORTS; ++sort)
		{
			for (int eta_bin=0;eta_bin<NetaBins;eta_bin++){
				sprintf(title,"p_{T}-efficiency after cuts for %s, %i ", sorts_of_particles[sort].Data(),centralityBin);
				h_pt_eff_after_final[centralityBin][sort][eta_bin] = new TH1F(Form("h_pt_eff_after_final%i%i%i",centralityBin,sort,eta_bin),title,100,0.,3.5);
				h_pt_eff_after_final[centralityBin][sort][eta_bin] -> Sumw2();
				f_pt_eff_after_final[centralityBin][sort][eta_bin] = new TF1(Form("f_pt_eff_after_final%i%i%i",centralityBin,sort,eta_bin),"pol8",0.2,2.);
				
				sprintf(title,"p_{T}^{reco} spectre after cuts for %s, %i ", sorts_of_particles[sort].Data(),centralityBin);
				h_pt_reco_final[centralityBin][sort][eta_bin] = new TH1F(Form("h_pt_reco_final%i%i%i",centralityBin,sort,eta_bin),title,100,0.,3.5);
				h_pt_reco_final[centralityBin][sort][eta_bin] -> Sumw2();
				
				sprintf(title,"p_{T}^{true} spectre after cuts for %s, %i ", sorts_of_particles[sort].Data(),centralityBin);
				h_pt_true_final[centralityBin][sort][eta_bin] = new TH1F(Form("h_pt_true_final%i%i%i",centralityBin,sort,eta_bin),title,100,0.,3.5);
				h_pt_true_final[centralityBin][sort][eta_bin] -> Sumw2();
			}
		}
	}

	for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin)
	{
		for (int sort = 0; sort < _N_SORTS; ++sort)
		{
			for (int eta_bin=0;eta_bin<NetaBins;eta_bin++){
				h_pt_reco_before[centralityBin][sort][eta_bin] = (TH1F*) inFile->Get(Form("h_pt[%i][%i][%i]",centralityBin,eta_bin,sort));
				h_pt_true_before[centralityBin][sort][eta_bin] = (TH1F*) inFile->Get(Form("h_pt_mc[%i][%i][%i]",centralityBin,eta_bin,sort));
				h_pt_reco_after[centralityBin][sort][eta_bin]  = (TH1F*) inFile->Get(Form("h_pt_after[%i][%i][%i]",centralityBin,eta_bin,sort));
				h_pt_true_after[centralityBin][sort][eta_bin]  = (TH1F*) inFile->Get(Form("h_pt_mc_after[%i][%i][%i]",centralityBin,eta_bin,sort));
				
				//sprintf(title,"p_{T}-efficiency before cuts for %s, %.2f - %.2f % ", sorts_of_particles[sort].Data(), centralityBinsRes[centralityBin], centralityBinsRes[centralityBin+1]);
				sprintf(title,"p_{T}-efficiency before cuts for %s, %i %i", sorts_of_particles[sort].Data(),centralityBin,eta_bin);

				//h_pt_eff_before[centralityBin][sort][eta_bin] = new TH1F(Form("h_pt_eff_before%i%i%i",centralityBin,sort,eta_bin),title,100,0.,3.5);
				//h_pt_eff_before[centralityBin][sort][eta_bin] -> Sumw2();
				//sprintf(title,"p_{T}-efficiency after cuts for %s, %2.0f-%2.0f % ", sorts_of_particles[sort].Data(), centralityBinsRes[centralityBin], centralityBinsRes[centralityBin+1]);
				sprintf(title,"p_{T}-efficiency after cuts for %s, %i %i", sorts_of_particles[sort].Data(), centralityBin,eta_bin);

				//h_pt_eff_after[centralityBin][sort][eta_bin]  = new TH1F(Form("h_pt_eff_after%i%i%i" ,centralityBin,sort,eta_bin),title,100,0.,3.5);
				//h_pt_eff_after[centralityBin][sort][eta_bin] -> Sumw2();
				
				//f_pt_eff_before[centralityBin][sort][eta_bin] = new TF1(Form("f_pt_eff_before%i%i%i",centralityBin,sort,eta_bin),"pol8",0.2,2.);
				//f_pt_eff_after[centralityBin][sort][eta_bin]  = new TF1(Form("f_pt_eff_before%i%i%i",centralityBin,sort,eta_bin),"pol8",0.2,2.);
			}
		}
	}*/

	for (int harm = 0; harm < _N_HARM; ++harm)
	{
		for (int _harm = 0; _harm < _N_HARM; ++_harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				sprintf(name,"p_Res2Psi_vs_b[%i][%i][%i]",harm,_harm,method);
				p_Res2Psi_vs_b[harm][_harm][method] = (TProfile*)inFile->Get(name);
				
				sprintf(name,"p_ResPsi_vs_b[%i][%i][%i]",harm,_harm,method);
				sprintf(title,"Res(#Psi_{%i,%s}) for v_{%i};b,fm;",_harm+1,methods_names[method].Data(),harm+1);
				p_Res[harm][_harm][method] = new TH1D(name,title,NcentralityBinsRes,centralityBinsRes);
			}
		}
	}
	
	outFile->cd();
	for (int harm = 0; harm < _N_HARM; ++harm)
	{
		for (int _harm = 0; _harm < _N_HARM; ++_harm)
		{
			for (int method = 0; method < _N_METHOD; ++method)
			{
				for (int centralitybin = 0; centralitybin < NcentralityBinsRes; ++centralitybin)
				{
                                    if (method == 0 && _harm == 1 && harm == 1)
                                    {
					Double_t res2 = p_Res2Psi_vs_b[harm][_harm][method]->GetBinContent(centralitybin + 1);
					if (res2 < 0)res2 = 0;
                                        p_Res[harm][_harm][method]->Fill(centralityBinsRes[centralitybin] + 0.1,TMath::Sqrt(res2));
                                    }
                                    else
                                    {
					Double_t res2 = p_Res2Psi_vs_b[harm][_harm][method]->GetBinContent(centralitybin + 1);
					if (res2 < 0)res2 = 0;
					Double_t chi_half = Chi(TMath::Sqrt(res2),harm+1);
					p_Res[harm][_harm][method]->Fill(centralityBinsRes[centralitybin] + 0.1,ResEventPlane(TMath::Sqrt(2.) * chi_half,harm+1));
					//p_Res[harm][_harm][method]->SetBinContent(centralitybin+1,ResEventPlane(TMath::Sqrt(2.) * chi_half,harm+1));
                                    }
				}
				
				char name[200];
				sprintf(name,"resolution_fit[%i][%i][%i]",harm,_harm,method);
				resolution_fit[harm][_harm][method] = new TF1(name,"pol8",0.01,100);
				Convert2Graph(p_Res[harm][_harm][method])->Fit(name,"W");
				p_Res[harm][_harm][method]->Write();
				resolution_fit[harm][_harm][method]->Write();
			}
		}
	}
	
	//for (Int_t eta_bin=0; eta_bin<NetaBins; eta_bin++){
		//for (int sort = 0; sort < _N_SORTS; ++sort){
			//h_pt_PID_eff[eta_bin][sort] -> Divide(h_pt_PID_efficiency_after_sorts[eta_bin][sort],h_pt_PID_efficiency_before[eta_bin],1.,1.);
			//h_pt_PID_eff[eta_bin][sort] -> Fit(Form("f_pt_PID_eff%i%i",eta_bin,sort),"WR");
			//h_pt_PID_eff[eta_bin][sort] -> Write();
			//f_pt_PID_eff[eta_bin][sort] -> Write();
		//}
	//}
	/*
	for (int centralityBin = 0; centralityBin < NcentralityBinsRes; ++centralityBin)
	{
		for (int sort = 0; sort < _N_SORTS; ++sort)
		{
			for (int eta_bin=0;eta_bin<NetaBins;eta_bin++){
				//h_pt_eff_after[centralityBin][sort][eta_bin] -> Divide(h_pt_reco_after[centralityBin][sort][eta_bin],h_pt_true_after[centralityBin][sort][eta_bin],1.,1.);
				//h_pt_eff_after[centralityBin][sort][eta_bin] -> Fit(Form("f_pt_eff_before%i%i%i",centralityBin,sort,eta_bin),"R");
				//h_pt_eff_after[centralityBin][sort] -> Write();
				//f_pt_eff_after[centralityBin][sort] -> Write();
				
				if (centralityBin >= 2 && centralityBin < 4) { 
					h_pt_reco_final[0][sort][eta_bin] ->Add(h_pt_reco_after[centralityBin][sort][eta_bin]);
					h_pt_true_final[0][sort][eta_bin] ->Add(h_pt_true_after[centralityBin][sort][eta_bin]);
				}
				if (centralityBin >= 4 && centralityBin < 8) { 
					h_pt_reco_final[1][sort][eta_bin] ->Add(h_pt_reco_after[centralityBin][sort][eta_bin]); 
					h_pt_true_final[1][sort][eta_bin] ->Add(h_pt_true_after[centralityBin][sort][eta_bin]);
				}
				if (centralityBin >= 8 && centralityBin < 10){ 
					h_pt_reco_final[2][sort][eta_bin] ->Add(h_pt_reco_after[centralityBin][sort][eta_bin]); 
					h_pt_true_final[2][sort][eta_bin] ->Add(h_pt_true_after[centralityBin][sort][eta_bin]);
				}
			}
		}
	}
	
	for (int centralityBin = 0; centralityBin < NcentralityBinsFlow; ++centralityBin)
	{
		for (int sort = 0; sort < _N_SORTS; ++sort)
		{
			for (int eta_bin=0;eta_bin<NetaBins;eta_bin++){
				h_pt_eff_after_final[centralityBin][sort][eta_bin] -> Divide(h_pt_reco_final[centralityBin][sort][eta_bin],h_pt_true_final[centralityBin][sort][eta_bin],1.,1.);
				h_pt_eff_after_final[centralityBin][sort][eta_bin] -> Fit(Form("f_pt_eff_after_final%i%i%i",centralityBin,sort,eta_bin),"R");
				h_pt_eff_after_final[centralityBin][sort][eta_bin] -> Write();
				f_pt_eff_after_final[centralityBin][sort][eta_bin] -> Write();
			}
		}
	}
*/
	outFile->Close();
	inFile->Close();
}
