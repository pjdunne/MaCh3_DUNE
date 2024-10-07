#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/samplePDFDUNEBeamFDBase.h"

int main(int argc, char * argv[]) {
  if(argc == 1){
    std::cout << "Usage: bin/EventRatesDUNEBeam config.cfg" << std::endl;
    return 1;
  }

  manager *fitMan = new manager(argv[1]);

  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs
  
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneBeamInstance(fitMan, DUNEPdfs, xsec, osc); 
 }
