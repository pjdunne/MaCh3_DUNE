#include "manager/manager.h"
#include "OscClass/OscClass_CUDAProb3.h"
#include "samplePDFDUNE/samplePDFDUNEBase.h"

int main(int argc, char * argv[]) {
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
  return 0;

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

  xsec->setParameters(XsecParVals);

  //####################################################################################
  
  std::vector<samplePDFDUNEBase*> SamplePDFs;
  for (std::string ConfigName : SampleConfigs) {
    samplePDFDUNEBase* SampleConfig = new samplePDFDUNEBase(POT, ConfigName.c_str(), xsec);
    
    SampleConfig->SetOscillator(Oscill);

    SamplePDFs.push_back(SampleConfig);
  }

  for (unsigned sample_i = 0 ; sample_i < SamplePDFs.size() ; ++sample_i) {
    SamplePDFs[sample_i] -> reweight(osc->getPropPars());
    std::cout << SamplePDFs[sample_i]->GetSampleName() << " : " << (SamplePDFs[sample_i]->get1DHist())->Integral() << std::endl;
  }

}
