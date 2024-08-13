#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/samplePDFDUNEBeamFDBase.h"
#include "samplePDFDUNE/samplePDFDUNEBeamNDBase.h"

int main(int argc, char * argv[]) {

  std::string ConfigName = argv[1];
  manager *FitManager = new manager(ConfigName.c_str());
  
  std::vector<std::string> xsecCovMatrixFile = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();

  std::cout << "==============================================================================" << std::endl;
  covarianceXsec *Xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");

  std::cout << "==============================================================================" << std::endl;
  std::string  OscMatrixFile = FitManager->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = FitManager->raw()["General"]["Systematics"]["OscCovName"].as<std::string>(); 
  std::vector<double> oscpars = FitManager->raw()["General"]["OscillationParameters"].as<std::vector<double>>();

  covarianceOsc* Osc = new covarianceOsc(OscMatrixName.c_str(), OscMatrixFile.c_str());

  std::cout << "==============================================================================" << std::endl;

  std::string OscillatorCfgName = FitManager->raw()["General"]["OscillatorConfigName"].as<std::string>();
  Oscillator* Oscill = new Oscillator(OscillatorCfgName);

  std::cout << "==============================================================================" << std::endl;

  samplePDFDUNEBeamFDBase *numu_FD_pdf = new samplePDFDUNEBeamFDBase(1.3628319e+23, "configs/AtmSample_numuselec.yaml", Xsec);
  samplePDFDUNEBeamNDBase *numu_ND_pdf = new samplePDFDUNEBeamNDBase(1.3628319e+23, "configs/AtmSample_numuselec.yaml", Xsec);
  numu_FD_pdf->SetOscillator(Oscill);
  
  std::cout << "==============================================================================" << std::endl;
  
  numu_FD_pdf->reweight(Osc->getPropPars());
  TH2D *numu_FD_asimov = (TH2D*)numu_FD_pdf->get2DHist()->Clone("numu_asimov");
  std::cout << numu_FD_asimov->Integral()/1. << std::endl;
}


