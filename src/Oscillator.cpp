#include "manager/manager.h"
#include "OscClass/OscClass_CUDAProb3.h"
#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include <fenv.h>

std::string getNameNoExt(std::string name, std::string ext)  
{                                                            
  std::size_t pos ;                                          
  pos = name.find(ext);                                      
  name = name.substr(0,pos);                                 
  return name ;                                              
}        

void saveCanvas(TCanvas* canvas, std::string name, std::string legend)                                                                  
{                                                            
  // name = getNameNoExt(name, ".root") ;                       
  // name = name + legend + ".root" ;                     
  // canvas -> SaveAs(name.c_str()) ;                           
                                                             
  // name = getNameNoExt(name, ".root") ;                       
  // name = name + ".png" ;                                     
  // canvas -> SaveAs(name.c_str()) ;                           
                                                             
  // name = getNameNoExt(name, ".png") ;                        
  // name = name + ".pdf" ;                                     
  // canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".pdf") ;                        
  name = name + ".eps" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
} 

int main(int argc, char * argv[]) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  manager *fitMan = new manager(argv[1]);
  
  std::string OscillatorCfgName = fitMan->raw()["General"]["OscillatorConfigName"].as<std::string>();
  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>();
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>();
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
  double POT = fitMan->raw()["General"]["POT"].as<double>();
  bool XsecParsAtGen = fitMan->raw()["General"]["Systematics"]["XsecAtGen"].as<bool>();
  std::vector<std::string> SampleConfigs = fitMan->raw()["General"]["SampleConfigs"].as<std::vector<std::string>>();

  //####################################################################################
  //Create Oscillation Covariance Matrix
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Creating Oscillation Covariance Matrix" << std::endl;
  covarianceOsc *osc = new covarianceOsc(OscMatrixName.c_str(),OscMatrixFile.c_str());

  // Taken from: https://github.com/DUNE/MaCh3_DUNE/blob/ce0ffaaf2919ef01f43a8a4103ac143e18345ad1/src/EventRatesDUNE.cpp#LL90C3-L90C88
  std::vector<double> oscpars{0.307,0.528,0.0218,7.53e-5, 2.509e-3,-1.601}; 
  osc->setParameters(oscpars);
  osc->acceptStep();

  //####################################################################################
  //Create Oscillator
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Creating Oscillator object" << std::endl;
  Oscillator* Oscill = new Oscillator(OscillatorCfgName);

  Oscill->FillOscillogram(osc->getPropPars(),25.0,0.5);
  Oscill->SaveOscillogramsToFile("oscillogram.root");

  //####################################################################################
  //Create XSec Covariance Matrix
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Creating XSec Covariance Matrix" << std::endl;
  covarianceXsec* xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str());

  // Taken from: https://github.com/DUNE/MaCh3_DUNE/blob/ce0ffaaf2919ef01f43a8a4103ac143e18345ad1/src/EventRatesDUNE.cpp#L144-L163
  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  if(XsecParsAtGen){
    TFile* XsecFile = new TFile(XsecMatrixFile.c_str(), "READ");
    TVectorD* XsecGeneratedParamArray = (TVectorD*)XsecFile->Get("xsec_param_nom");
    std::cout << "Setting xsec systs to their generated values " << std::endl;
    for (unsigned param_i = 0 ; param_i < XsecParVals.size() ; ++param_i) {
      std::cout << "Generated value for param " << param_i << " is " << (*XsecGeneratedParamArray)(param_i) << std::endl;
      XsecParVals[param_i] = (*XsecGeneratedParamArray)(param_i);
      std::cout << "Set parameter " << param_i << " to value " << XsecParVals[param_i] << std::endl;
    }
  }
  else{
    std::cout << "Keeping xsec parameters at their prior values" << std::endl;
  }

  // 
  xsec->setParameters(XsecParVals);

  //####################################################################################
  std::vector<samplePDFDUNEBase*> SamplePDFs;
  for (std::string ConfigName : SampleConfigs) {
    samplePDFDUNEBase* SampleConfig = new samplePDFDUNEBase(POT, ConfigName.c_str(), xsec);
    
    SampleConfig->SetOscillator(Oscill);

    SamplePDFs.push_back(SampleConfig);
  }

  // unosc
  std::vector<double> oscpars_un(oscpars);
  oscpars_un[0] = 0;
  oscpars_un[1] = 0;
  oscpars_un[2] = 0;
  oscpars_un[3] = 0;
  oscpars_un[4] = 0;

//Some place to store the histograms
  std::vector<TH2D*> oscillated_hists;
  std::vector<TH2D*> unoscillated_hists;
  std::vector<std::string> sample_names;

  TFile* Outfile = new TFile("EventsRateAtm.root", "RECREATE");

  for (unsigned sample_i = 0 ; sample_i < SamplePDFs.size() ; ++sample_i) {
    
	
	std::string name = SamplePDFs[sample_i]->GetSampleName();
	sample_names.push_back(name);
	TString NameTString = TString(name.c_str());
	// Unoscillated
	osc -> setParameters(oscpars_un);
	osc -> acceptStep();

	SamplePDFs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
	SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	TH2D *numu_unosc = (TH2D*)SamplePDFs[sample_i] -> get2DHist() -> Clone(NameTString+"_unosc");
	unoscillated_hists.push_back(numu_unosc);

	TCanvas *nomcanv = new TCanvas("nomcanv","",1200,600);
	numu_unosc -> Draw("HIST");

	std::string plotname = "Dummy_Hist" ;
	saveCanvas(nomcanv, plotname,"_nominal_spectra") ;
	Outfile->cd();
	numu_unosc->Write(NameTString);

	osc->setParameters(oscpars);
	osc->acceptStep();
	// Oscillated
	std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	  << "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	  << "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	  << "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	  << "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	  << "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;

	SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	TH2D *numu_osc = (TH2D*)SamplePDFs[sample_i] -> get2DHist()->Clone(NameTString+"_osc");
	oscillated_hists.push_back(numu_osc);

  }

  //Now print out some event rates, we'll make a nice latex table at some point 
  for (unsigned sample_i = 0; sample_i < SamplePDFs.size() ; ++sample_i) {
	std::cout << "Integrals of nominal hists: " << std::endl;
	std::cout << sample_names[sample_i].c_str() << " unosc:      " << unoscillated_hists[sample_i]-> Integral() << std::endl;
	std::cout << sample_names[sample_i].c_str() << "   osc:      " << oscillated_hists[sample_i]-> Integral() << std::endl; 
	std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
  }

  Outfile->Write();
  Outfile->Close();




  // for (unsigned sample_i = 0 ; sample_i < SamplePDFs.size() ; ++sample_i) {
  //   SamplePDFs[sample_i] -> reweight(osc->getPropPars());
  //   std::cout << SamplePDFs[sample_i]->GetSampleName() << " : " << (SamplePDFs[sample_i]->get1DHist())->Integral() << std::endl;
  //   SamplePDFs[sample_i]->get1DHist()->Dump();
    
  //   TH1D *numu_unosc = (TH1D*)SamplePDFs[sample_i] -> get1DHist() -> Clone("unosc");
  //   std::cout << "nb entries: " << numu_unosc->GetEntries() << std::endl;

  //   TCanvas *nomcanv = new TCanvas("nomcanv","",1200,600);
  //   numu_unosc -> Draw("HIST");

  //   std::string plotname = "Dummy_Hist_" ;
  //   saveCanvas(nomcanv, plotname,"_nominal_spectra") ;
  // }

  // // Oscillated
	// std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	//   << "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	//   << "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	//   << "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	//   << "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	//   << "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;


	// SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	// TH1D *numu_osc = (TH1D*)SamplePDFs[sample_i] -> get1DHist()->Clone(NameTString+"_osc");
	// oscillated_hists.push_back(numu_osc);

  // }

  // //Now print out some event rates, we'll make a nice latex table at some point 
  // for (unsigned sample_i = 0; sample_i < SamplePDFs.size() ; ++sample_i) {
	// std::cout << "Integrals of nominal hists: " << std::endl;
	// std::cout << sample_names[sample_i].c_str() << " unosc:      " << unoscillated_hists[sample_i]-> Integral() << std::endl;
	// std::cout << sample_names[sample_i].c_str() << "   osc:      " << oscillated_hists[sample_i]-> Integral() << std::endl; 
	// std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
  // }

  return 0;

}
