#include "TROOT.h"
#include "TString.h"
#include <iostream>

#include "MpdCalculator.h"

R__ADD_LIBRARY_PATH("/lustre/nyx/hades/user/parfenov/mpd_new/real-flow") // if needed
R__LOAD_LIBRARY(MpdCalculator_C.so)
R__LOAD_LIBRARY(utility_C.so)

int main_flow(TString inFileName , TString outFileName, TString resFitFile, TString dcaFile)
{
	//gSystem->Load("./utility_C.so");
	//gSystem->Load("./MpdCalculator_C.so");
	
	MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFile);
	mpd.CalculateFlow(0, resFitFile.Data());
	mpd.Write();
        return 0;
}
