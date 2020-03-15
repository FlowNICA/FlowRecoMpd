//#include "MpdCalculator.C"
//#include "utility.C"

#include "TROOT.h"
#include "TString.h"
#include <iostream>

#include "MpdCalculator.h"

R__ADD_LIBRARY_PATH("/lustre/nyx/hades/user/parfenov/mpd_new/real-flow") // if needed
R__LOAD_LIBRARY(MpdCalculator_C.so)
R__LOAD_LIBRARY(utility_C.so)


int main_resolutions(TString inFileName , TString outFileName, TString dcaFileName)
{
  std::cout << "loading libraries" << std::endl;
  //gSystem->Load("./utility_C.so");
  //gSystem->Load("/lustre/nyx/hades/user/parfenov/real-flow/utility_C.so");
  //gSystem->Load("./MpdCalculator_C.so");

  MpdCalculator mpd = MpdCalculator(inFileName,outFileName,dcaFileName);
  mpd.CalculateResolutions(0);
  mpd.Write();
  return 0;
}
