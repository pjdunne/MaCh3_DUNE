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

#include "mcmc/mcmc.h"
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"

int main(int argc, char * argv[]) {
  if(argc == 1){
    MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
    return 1;
  }
  manager* FitManager = new manager(argv[1]);

  //###############################################################################################################################
  //Create samplePDFFD objects
  
  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;
  
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);

  //###############################################################################################################################
  //Perform reweight, print total integral, and set data

  std::vector<TH1*> DUNEHists;
  for(auto Sample : DUNEPdfs){
    Sample->reweight();
    if (Sample->GetNDim() == 1)
      DUNEHists.push_back(Sample->get1DHist());
    else if (Sample->GetNDim() == 2)
      DUNEHists.push_back(Sample->get2DHist());

    MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetName(), Sample->get1DHist()->Integral());
    if (Sample->GetNDim() == 1)
      Sample->addData((TH1D*)DUNEHists.back());
    else if (Sample->GetNDim() == 2)
      Sample->addData((TH2D*)DUNEHists.back());
  }
  std::unique_ptr<FitterBase> MaCh3Fitter = std::make_unique<mcmc>(FitManager);

  //###############################################################################################################################
  //Lets benefit from the core code utilities 
  
  //Add samples to FitterBase
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->addSamplePDF(Sample);
  }

  MaCh3Fitter->addSystObj(osc);
  MaCh3Fitter->addSystObj(xsec);
  
  MaCh3Fitter->RunLLHScan();  
}
