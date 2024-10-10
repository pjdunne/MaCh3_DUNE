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

void Write1DHistogramsToFile(manager* &Manager, std::vector<TH1D*> Histograms) {

}

int main(int argc, char * argv[]) {
  if(argc == 1){
    std::cout << "Usage: bin/EventRatesDUNEBeam config.cfg" << std::endl;
    return 1;
  }

  manager *fitMan = new manager(argv[1]);

  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFFD objects
  
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec, osc); 


  std::vector<TH1D*> DUNEHists;
  for(auto Sample : DUNEPdfs){
	Sample->reweight();
    DUNEHists.push_back(Sample->get1DHist());

    std::string EventRateString = fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetName(), EventRateString);

  }

  // Now write out the saved hsitograms to file 
  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");

  auto OutputFile = std::unique_ptr<TFile>(new TFile(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for(auto Hists : DUNEHists){
    Hists->Write();
  }
  OutputFile->Close();
}
