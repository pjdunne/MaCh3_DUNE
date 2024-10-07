#include "samplePDFDUNE/MaCh3DUNEFactory.h"

#include "samplePDFDUNE/samplePDFDUNEBeamFDBase.h"

void MakeMaCh3DuneBeamInstance(manager *FitManager, std::vector<samplePDFFDBase*> &DUNEPdfs, covarianceXsec *&xsec, covarianceOsc *&osc){

  // there's a check inside the manager class that does this; left here for demonstrative purposes
  if (FitManager == nullptr) {
    std::cerr << "Didn't find a good config in input configuration" << std::endl;
    throw;
  }
  
  //Check that you have specified some DUNE samples
  if(!FitManager->raw()["General"]["DUNESamples"]){
    std::cerr << "You didn't specify any DUNESample configs to create samples from. Please add General:DUNESamples to your config" << std::endl;
    throw;
  }
  
  //Check that you have specified some DUNE samples
  if(!FitManager->raw()["General"]["DUNESamplesPOT"]){
    std::cerr << "You didn't specify the POT for each DUNE sample configs to create samples from. Please add General:DUNESamplesPOT to your config" << std::endl;
    throw;
  }
  
  // Get inputted systematic parameters covariance matrices
  std::vector<std::string> xsecCovMatrixFile = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();
  
  // Setup the covariance matrices
  if(xsec == nullptr){
    xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
  }
  else{
    std::cout << "covariance Xsec has already been created so I am not re-initialising the object" << std::endl; 
  }
  std::cout << "cov xsec setup" << std::endl;
  std::cout << "------------------------------" << std::endl;
  
  std::string  OscMatrixFile = FitManager->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = FitManager->raw()["General"]["Systematics"]["OscCovName"].as<std::string>(); 
  std::vector<double> oscpars = FitManager->raw()["General"]["OscillationParameters"].as<std::vector<double>>();
  std::cout << "Oscillation Parameters: { ";
  for (unsigned int i=0;i<oscpars.size();i++) {
    std::cout << oscpars[i] << " ";
  }
  std::cout << "}" << std::endl;
  
  osc = new covarianceOsc(OscMatrixName.c_str(), OscMatrixFile.c_str());
  osc->setName("osc_cov");
  std::cout << "Osc cov setup " << std::endl;
  std::cout << "------------------------------" << std::endl;
  
  // ==========================================================
  //read flat prior, fixed paramas from the config file
  std::vector<std::string> XsecFixParams   = GetFromManager<std::vector<std::string>>(FitManager->raw()["General"]["Systematics"]["XsecFix"], {""});
  
  // Fixed xsec parameters loop
  if (XsecFixParams.size() == 1 && XsecFixParams.at(0) == "All") {
    for (int j = 0; j < xsec->getSize(); j++) {
      xsec->toggleFixParameter(j);
    }
  } else {
    for (unsigned int j = 0; j < XsecFixParams.size(); j++) {
      xsec->toggleFixParameter(XsecFixParams.at(j));
    }
  }
  std::cout << "xsec parameters loop done" << std::endl;
  
  // Fill the parameter values with their nominal values
  // should _ALWAYS_ be done before overriding with fix or flat
  xsec->setParameters();
  osc->setParameters(oscpars); 
  
  //####################################################################################
  //Create samplePDFDUNEBase Objs
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "Loading T2K samples.." << "\n" << std::endl;
  std::vector<std::string> DUNESampleConfigs = FitManager->raw()["General"]["DUNESamples"].as<std::vector<std::string>>();
  std::vector<double> DUNESamplePOTs = FitManager->raw()["General"]["DUNESamplesPOT"].as<std::vector<double>>();
  
  if(DUNESampleConfigs.size() != DUNESamplePOTs.size()){
    std::cerr << "Size of DUNESamples and DUNESamplesPOT is not the same and they need to be!" << std::endl;
    throw;
  }
  
  for(unsigned int Sample_i = 0 ; Sample_i < DUNESampleConfigs.size() ; Sample_i++){

    manager* tempSampleManager = new manager(DUNESampleConfigs[Sample_i].c_str());
    std::string SampleType = tempSampleManager->raw()["SampleType"].as<std::string>();

    samplePDFFDBase* FDBaseSample;
    if (SampleType == "BeamFD") {
      samplePDFDUNEBeamFDBase* BeamFDSample = new samplePDFDUNEBeamFDBase(DUNESamplePOTs[Sample_i], DUNESampleConfigs[Sample_i].c_str(), xsec);
      FDBaseSample = (samplePDFFDBase*)BeamFDSample;
    } else {
      std::cerr << "Invalid SampleType:" << SampleType << " defined in config:" << DUNESampleConfigs[Sample_i] << std::endl;
      throw;
    }

    FDBaseSample->UseNonDoubledAngles(true);
    FDBaseSample->SetXsecCov(xsec);
    FDBaseSample->SetOscCov(osc);
    FDBaseSample->SetupOscCalc(osc->GetPathLength(), osc->GetDensity());

    // Pure for debugging, lets us set which weights we don't want via the manager
#if DEBUG_DUNE_WEIGHTS==1
    FDBaseSample->setWeightSwitchOffVector(FitManager->getWeightSwitchOffVector());
    FDBaseSample->setXsecWeightSwitchOffVector(FitManager->getXsecWeightSwitchOffVector());
#endif

    DUNEPdfs.push_back(FDBaseSample);
  }
  
  return;
}
