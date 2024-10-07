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
  // ----------------------- OPTIONS ---------------------------------------- //
  if(argc == 1){
    std::cout << "Usage: bin/EventRatesDUNEBeam config.cfg" << std::endl;
    return 1;
  }

  manager *fitMan = new manager(argv[1]);
  auto OutputFileName = fitMan->raw()["General"]["OutputFile"].as<std::string>();

  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs
  std::cout << "Loading DUNE samples.." << "\n" << std::endl;
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneBeamInstance(fitMan, DUNEPdfs, xsec, osc); 
  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  xsec->setParameters(XsecParVals);
  xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());

  //Some place to store the histograms
  std::vector<TH1D*> oscillated_hists;
  std::vector<TH1D*> unoscillated_hists;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(new TFile(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {
    std::string name = DUNEPdfs[sample_i]->GetName();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    
    DUNEPdfs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
    osc->setParameters();
    DUNEPdfs[sample_i] -> reweight();
    TH1D *SampleHistogram = (TH1D*)DUNEPdfs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc");
    unoscillated_hists.push_back(SampleHistogram);
  }
 
  //Now print out some event rates, we'll make a nice latex table at some point 
  for (unsigned iPDF = 0; iPDF < DUNEPdfs.size() ; ++iPDF) {
	std::cout << "Integrals of nominal hists: " << std::endl;
	std::cout << sample_names[iPDF].c_str() << ": " << unoscillated_hists[iPDF]-> Integral() << std::endl;
	std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
  }

  OutputFile->Write();
  OutputFile->Close();
  return 0;
 }
